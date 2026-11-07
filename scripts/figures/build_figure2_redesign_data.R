if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(digest)
})

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

root <- Sys.getenv("MHFNDD_PROJECT_ROOT")
out_dir <- file.path(root, "FIGURE2_REDESIGN_DATA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

src <- c(
  frozen_loadings = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/BIDIRECTIONAL_SPECIES_LOADING_H2.tsv"),
  rbo = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/GLOBAL_RBO_ALL_PAIRS.tsv"),
  topk = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/GLOBAL_TOPK_OVERLAP_CURVES.tsv"),
  core_species = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/GLOBAL_COMPOSITION_CORE_SPECIES.tsv"),
  permutation_script = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/scripts/01_phase7_global_composition.R"),
  null_calibration = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/GLOBAL_RBO_NULL_CALIBRATION.md"),
  phase7_hashes = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/logs/formal_output_sha256.tsv"),
  taxonomy = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE8_PHYLOGENETIC_ORGANIZATION/PHYLOGENETIC_SPECIES_MAPPING.tsv"),
  final_workbook = file.path(root, "SUBMISSION_PACKAGE/SUPPLEMENTARY_TABLES_FINAL/Supplementary_Tables.xlsx"),
  supplement_crosswalk = file.path(root, "SUBMISSION_PACKAGE/SUPPLEMENTARY_TABLES_FINAL/SUPPLEMENTARY_TABLE_SOURCE_CROSSWALK.tsv"),
  v3_manuscript = file.path(root, "MANUSCRIPT_INTEGRATION_FREEZE_20261031/REVISION_3/GLOBAL_GFA_INTEGRATED_DRAFT_V3.md"),
  v3_handoff = file.path(root, "MANUSCRIPT_INTEGRATION_FREEZE_20261031/REVISION_3/MASTER_RESEARCH_HANDOFF_V3.md"),
  evidence_ledger = file.path(root, "MANUSCRIPT_INTEGRATION_FREEZE_20261031/REVISION_2/EVIDENCE_LEDGER_V2.tsv")
)
stopifnot(all(file.exists(src)))

read_char_tsv <- function(path) {
  read.delim(path, colClasses = "character", check.names = FALSE, quote = "", comment.char = "", na.strings = c("NA", "<NA>"))
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, col.names = TRUE,
              quote = FALSE, na = "")
}

num <- function(x) suppressWarnings(as.numeric(x))
bool_chr <- function(x) ifelse(x, "TRUE", "FALSE")
direction_chr <- function(x) {
  z <- num(x)
  ifelse(z > 0, "POSITIVE", ifelse(z < 0, "NEGATIVE", "ZERO"))
}

key_string <- function(df, keys) {
  apply(as.data.frame(df[, keys, drop = FALSE]), 1, function(z) paste(ifelse(is.na(z), "<NA>", z), collapse = "|"))
}

compare_frames <- function(observed, expected, keys, tolerance = 1e-10) {
  observed <- as.data.frame(observed, stringsAsFactors = FALSE)
  expected <- as.data.frame(expected, stringsAsFactors = FALSE)
  stopifnot(all(names(expected) %in% names(observed)))
  observed <- observed[, names(expected), drop = FALSE]
  ko <- key_string(observed, keys)
  ke <- key_string(expected, keys)
  stopifnot(length(ko) == length(ke), !anyDuplicated(ko), !anyDuplicated(ke), setequal(ko, ke))
  observed <- observed[match(ke, ko), , drop = FALSE]
  numeric_diffs <- numeric()
  character_mismatches <- integer()
  for (v in names(expected)) {
    if (is.numeric(expected[[v]])) {
      xo <- suppressWarnings(as.numeric(observed[[v]]))
      xe <- expected[[v]]
      stopifnot(identical(is.na(xo), is.na(xe)))
      if (any(!is.na(xe))) numeric_diffs[v] <- max(abs(xo[!is.na(xe)] - xe[!is.na(xe)]))
    } else {
      xo <- as.character(observed[[v]]); xe <- as.character(expected[[v]])
      xo[is.na(xo)] <- ""; xe[is.na(xe)] <- ""
      character_mismatches[v] <- sum(xo != xe)
    }
  }
  max_numeric_diff <- if (length(numeric_diffs)) max(numeric_diffs) else 0
  total_character_mismatches <- sum(character_mismatches)
  stopifnot(max_numeric_diff <= tolerance, total_character_mismatches == 0L)
  list(n = nrow(expected), max_numeric_diff = max_numeric_diff,
       character_mismatches = total_character_mismatches)
}

load_chr <- read_char_tsv(src[["frozen_loadings"]])
load_num <- read.delim(src[["frozen_loadings"]], check.names = FALSE)
tax_chr <- read_char_tsv(src[["taxonomy"]])
rbo_chr <- read_char_tsv(src[["rbo"]])
rbo_num <- read.delim(src[["rbo"]], check.names = FALSE)
topk_chr <- read_char_tsv(src[["topk"]])
topk_num <- read.delim(src[["topk"]], check.names = FALSE)
core_chr <- read_char_tsv(src[["core_species"]])
core_num <- read.delim(src[["core_species"]], check.names = FALSE)

required_worksheets <- c(
  "S3_Loadings_PVE", "S4A_Factor_Correlations", "S4B_Factor_Matching",
  "S5A_TopK_Overlap", "S5B_Shared_Species", "S6A_RBO_Permutation",
  "S6B_TopK_Curves"
)
workbook_sheet_names <- excel_sheets(src[["final_workbook"]])
stopifnot(all(required_worksheets %in% workbook_sheet_names))
wb <- setNames(lapply(required_worksheets, function(s) {
  read_excel(src[["final_workbook"]], sheet = s, skip = 4)
}), required_worksheets)

# Inputs listed in the manifest are content-read here; roles below describe the exact checks.
phase7_hash_manifest <- read.delim(src[["phase7_hashes"]], check.names = FALSE)
supplement_crosswalk <- read.delim(src[["supplement_crosswalk"]], check.names = FALSE)
permutation_script_text <- paste(readLines(src[["permutation_script"]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
null_calibration_text <- paste(readLines(src[["null_calibration"]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
v3_manuscript_text <- paste(readLines(src[["v3_manuscript"]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
v3_handoff_text <- paste(readLines(src[["v3_handoff"]], warn = FALSE, encoding = "UTF-8"), collapse = "\n")
evidence_ledger <- read.delim(src[["evidence_ledger"]], check.names = FALSE)

stopifnot(
  grepl("SEED <- 20261021L", permutation_script_text, fixed = TRUE),
  grepl("B <- 20000L", permutation_script_text, fixed = TRUE),
  grepl("P_GRID <- c(0.90, 0.95, 0.975)", permutation_script_text, fixed = TRUE),
  grepl("K_GRID <- c(5L, 10L, 15L, 20L, 25L, 30L, 40L, 50L, 60L, 80L)", permutation_script_text, fixed = TRUE),
  grepl("20000 deterministic permutations (seed 20261021)", null_calibration_text, fixed = TRUE),
  grepl("seed 20261022", null_calibration_text, fixed = TRUE),
  grepl("# Master research handoff V3", v3_handoff_text, fixed = TRUE),
  grepl("Revision 3", v3_handoff_text, fixed = TRUE)
)

# Verify that the frozen Phase7 hash ledger is not merely listed but agrees with
# the three numerical Phase7 files consumed below.
phase7_hash_targets <- c(
  GLOBAL_RBO_ALL_PAIRS.tsv = src[["rbo"]],
  GLOBAL_TOPK_OVERLAP_CURVES.tsv = src[["topk"]],
  GLOBAL_COMPOSITION_CORE_SPECIES.tsv = src[["core_species"]]
)
phase7_hash_rows <- phase7_hash_manifest[match(names(phase7_hash_targets), phase7_hash_manifest$file), ]
stopifnot(!anyNA(phase7_hash_rows$file))
phase7_actual_hashes <- vapply(phase7_hash_targets, function(x) {
  digest(file = x, algo = "sha256", serialize = FALSE)
}, character(1))
stopifnot(identical(unname(as.character(phase7_hash_rows$sha256)), unname(phase7_actual_hashes)))

# Verify that the final crosswalk explicitly assigns the seven consumed sheets
# to the expected frozen source families.
crosswalk_expected <- data.frame(
  worksheet = required_worksheets,
  supplementary_table = c("S3", "S4", "S4", "S5", "S5", "S6", "S6"),
  stringsAsFactors = FALSE
)
crosswalk_observed <- supplement_crosswalk[
  match(crosswalk_expected$worksheet, supplement_crosswalk$worksheet), ]
stopifnot(!anyNA(crosswalk_observed$worksheet),
          identical(as.character(crosswalk_observed$supplementary_table), crosswalk_expected$supplementary_table),
          all(grepl("Frozen", crosswalk_observed$controlling_status, fixed = TRUE)))

# These manuscript/ledger checks are deliberately limited to the rounded claims
# that the files actually contain. The unrounded numerical checks use TSV/XLSX.
manuscript_claim_checks <- c(
  grepl("RBO was 0.317", v3_manuscript_text, fixed = TRUE) &&
    grepl("0.303 for F2", v3_manuscript_text, fixed = TRUE),
  grepl("shared eight and seven top-20 species", v3_manuscript_text, fixed = TRUE),
  grepl("full 241-species loading correlations were only 0.231 and 0.419", v3_manuscript_text, fixed = TRUE),
  grepl("principal angles of 59.0", v3_manuscript_text, fixed = TRUE) &&
    grepl("and 80.9", v3_manuscript_text, fixed = TRUE)
)
stopifnot(all(manuscript_claim_checks))
ledger_expected <- c(
  E004 = "0.3168;0.3032", E005 = "0.00015;0.00085",
  E006 = "0.0349;0.1079", E007 = "8;7", E008 = "59.0;80.9"
)
ledger_observed <- setNames(as.character(evidence_ledger$numeric_value), evidence_ledger$evidence_id)
stopifnot(identical(unname(ledger_observed[names(ledger_expected)]), unname(ledger_expected)))

stopifnot(nrow(load_chr) == 241L, length(unique(load_chr$trait_order)) == 241L,
          length(unique(load_chr$species)) == 241L)

# Verify all 241 species and all four loading profiles against final S3.
s3 <- wb[["S3_Loadings_PVE"]]
names(s3) <- c(
  "trait_order", "swedish_accession", "hunt_accession", "species",
  "Swedish_F1_loading_raw", "Swedish_F1_abs_rank", "Swedish_F1_PVE",
  "Swedish_F2_loading_raw", "Swedish_F2_abs_rank", "Swedish_F2_PVE",
  "HUNT_F1_loading_raw", "HUNT_F1_abs_rank", "HUNT_F1_PVE",
  "HUNT_F2_loading_raw", "HUNT_F2_abs_rank", "HUNT_F2_PVE"
)
s3 <- s3 %>%
  mutate(trait_order = suppressWarnings(as.numeric(trait_order))) %>%
  filter(is.finite(trait_order)) %>%
  arrange(trait_order)
load_num <- load_num %>% arrange(trait_order)
stopifnot(nrow(s3) == 241L, identical(as.character(s3$species), load_num$species))
loading_fields <- c("Swedish_F1_loading_raw", "Swedish_F2_loading_raw", "HUNT_F1_loading_raw", "HUNT_F2_loading_raw")
rank_fields <- c("Swedish_F1_abs_rank", "Swedish_F2_abs_rank", "HUNT_F1_abs_rank", "HUNT_F2_abs_rank")
workbook_loading_max_abs_diff <- max(vapply(loading_fields, function(v) max(abs(as.numeric(s3[[v]]) - load_num[[v]])), numeric(1)))
workbook_rank_max_abs_diff <- max(vapply(rank_fields, function(v) max(abs(as.numeric(s3[[v]]) - load_num[[v]])), numeric(1)))
workbook_display_precision_only <- workbook_loading_max_abs_diff < 1e-10 && workbook_rank_max_abs_diff == 0
stopifnot(workbook_display_precision_only)

# Read the frozen HUNT sign multipliers from final S4B; do not define them independently.
s4b <- wb[["S4B_Factor_Matching"]]
names(s4b) <- c(
  "swedish_factor", "hunt_factor", "sign_applied_to_hunt",
  "species_pearson", "species_spearman", "snp_pearson", "snp_spearman",
  "species_pearson_aligned", "species_spearman_aligned",
  "snp_pearson_aligned", "snp_spearman_aligned"
)
s4b <- s4b %>%
  mutate(swedish_factor_numeric = suppressWarnings(as.numeric(swedish_factor))) %>%
  filter(is.finite(swedish_factor_numeric)) %>%
  mutate(swedish_factor = paste0("F", swedish_factor_numeric),
         hunt_factor = paste0("F", as.numeric(hunt_factor)),
         sign_applied_to_hunt = as.numeric(sign_applied_to_hunt))
stopifnot(nrow(s4b) == 2L, setequal(s4b$hunt_factor, c("F1", "F2")),
          !anyDuplicated(s4b$hunt_factor))
hunt_sign_map <- setNames(s4b$sign_applied_to_hunt, s4b$hunt_factor)
stopifnot(identical(unname(hunt_sign_map[c("F1", "F2")]), c(-1, 1)))

f2a_alignment_max_abs_diff <- max(c(
  abs(load_num$HUNT_F1_loading_raw * hunt_sign_map[["F1"]] - load_num$HUNT_loading_aligned_to_Swedish_F1),
  abs(load_num$HUNT_F2_loading_raw * hunt_sign_map[["F2"]] - load_num$HUNT_loading_aligned_to_Swedish_F2)
))
stopifnot(f2a_alignment_max_abs_diff < 1e-15,
          all(load_num$HUNT_sign_applied_for_Swedish_F1 == hunt_sign_map[["F1"]]),
          all(load_num$HUNT_sign_applied_for_Swedish_F2 == hunt_sign_map[["F2"]]))

# Verify Phase8 taxonomy rows retain the same Swedish raw loadings and frozen HUNT sign alignment.
tax_num <- read.delim(src[["taxonomy"]], check.names = FALSE)
tax_num <- tax_num %>% arrange(trait_order)
stopifnot(nrow(tax_num) == 241L, identical(tax_num$frozen_species, load_num$species))
phase8_loading_max_abs_diff <- max(c(
  abs(tax_num$Swedish_F1_signed_loading - load_num$Swedish_F1_loading_raw),
  abs(tax_num$Swedish_F2_signed_loading - load_num$Swedish_F2_loading_raw),
  abs(tax_num$HUNT_F1_signed_loading - load_num$HUNT_loading_aligned_to_Swedish_F1),
  abs(tax_num$HUNT_F2_signed_loading - load_num$HUNT_loading_aligned_to_Swedish_F2)
), na.rm = TRUE)
stopifnot(phase8_loading_max_abs_diff < 1e-12)

# Panel A: preserve source text for all original loading/rank values.
tax_keep <- tax_chr %>%
  select(trait_order, frozen_species, gtdb_tree_tip, frozen_gtdb_release,
         swedish_assembly_accession, hunt_hMGS, r232_manifest_formatted_species,
         mapping_class, mapping_basis, uniquely_mapped, domain, phylum, class,
         order, family, genus, species_taxon, reference_release)

f2a <- load_chr %>%
  select(
    trait_order, harmonized_species_name = species, swedish_accession, hunt_accession,
    Swedish_F1_loading_raw, Swedish_F1_abs_rank,
    Swedish_F2_loading_raw, Swedish_F2_abs_rank,
    HUNT_F1_loading_raw, HUNT_F1_abs_rank,
    HUNT_F2_loading_raw, HUNT_F2_abs_rank,
    HUNT_loading_aligned_to_Swedish_F1, HUNT_source_factor_for_Swedish_F1,
    HUNT_sign_applied_for_Swedish_F1,
    HUNT_loading_aligned_to_Swedish_F2, HUNT_source_factor_for_Swedish_F2,
    HUNT_sign_applied_for_Swedish_F2,
    Swedish_F1_top20, Swedish_F1_top40, Swedish_F2_top20, Swedish_F2_top40,
    HUNT_F1_top20, HUNT_F1_top40, HUNT_F2_top20, HUNT_F2_top40
  ) %>%
  mutate(
    Swedish_F1_direction_raw = direction_chr(Swedish_F1_loading_raw),
    Swedish_F2_direction_raw = direction_chr(Swedish_F2_loading_raw),
    HUNT_F1_direction_raw = direction_chr(HUNT_F1_loading_raw),
    HUNT_F2_direction_raw = direction_chr(HUNT_F2_loading_raw),
    HUNT_F1_direction_after_frozen_alignment = direction_chr(HUNT_loading_aligned_to_Swedish_F1),
    HUNT_F2_direction_after_frozen_alignment = direction_chr(HUNT_loading_aligned_to_Swedish_F2)
  ) %>%
  left_join(tax_keep, by = c("trait_order", "harmonized_species_name" = "frozen_species")) %>%
  arrange(num(trait_order))
stopifnot(nrow(f2a) == 241L)
write_tsv(f2a, file.path(out_dir, "F2A_factor_loadings.tsv"))

# Panel B: exact frozen cross-cohort RBO rows for all four factor pairs and all saved p values.
f2b <- rbo_chr %>%
  filter(comparison_scope == "cross_cohort_pair",
         pair_id %in% c("F1_F1", "F1_F2", "F2_F1", "F2_F2")) %>%
  mutate(
    primary_matrix_p = ifelse(rbo_p == "0.95", "TRUE", "FALSE"),
    permutation_count = "20000",
    primary_permutation_seed = "20261021",
    propensity_sensitivity_seed = "20261022",
    primary_null_definition = "HUNT row block permuted relative to Swedish identities within 9 strata formed by PC1/PC2 tertiles from 8 standardized information covariates",
    primary_null_stratum_sizes = "26;23;31;29;26;25;25;31;25",
    propensity_null_definition = "HUNT row block permuted within information-PC1 tertile by pre-existing non-global high-loading propensity class",
    propensity_null_stratum_sizes = "35;14;31;54;10;16;61;10;10",
    saved_rbo_p_grid = "0.90;0.95;0.975",
    saved_topk_grid = "5;10;15;20;25;30;40;50;60;80",
    source_phase = "PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY"
  )
stopifnot(nrow(f2b) == 12L)
write_tsv(f2b, file.path(out_dir, "F2B_RBO_matrix.tsv"))

# Helpers for all four Swedish-HUNT factor comparisons.
pair_defs <- data.frame(
  pair_id = c("F1_F1", "F1_F2", "F2_F1", "F2_F2"),
  swedish_factor = c("F1", "F1", "F2", "F2"),
  hunt_factor = c("F1", "F2", "F1", "F2"),
  hunt_raw_sign_multiplier = as.character(unname(hunt_sign_map[c("F1", "F2", "F1", "F2")])),
  pair_status = c("primary_frozen", "off_diagonal_sensitivity", "off_diagonal_sensitivity", "primary_frozen"),
  pairing_scheme = c("identity_frozen", "swapped_sensitivity", "swapped_sensitivity", "identity_frozen"),
  stringsAsFactors = FALSE
)

get_col <- function(df, prefix, factor, suffix) df[[paste0(prefix, "_", factor, "_", suffix)]]

rank_rows <- lapply(seq_len(nrow(pair_defs)), function(i) {
  p <- pair_defs[i, ]
  sw_load <- get_col(load_chr, "Swedish", p$swedish_factor, "loading_raw")
  sw_rank <- get_col(load_chr, "Swedish", p$swedish_factor, "abs_rank")
  hu_load <- get_col(load_chr, "HUNT", p$hunt_factor, "loading_raw")
  hu_rank <- get_col(load_chr, "HUNT", p$hunt_factor, "abs_rank")
  mult <- as.numeric(p$hunt_raw_sign_multiplier)
  data.frame(
    pair_id = p$pair_id,
    pair_status = p$pair_status,
    pairing_scheme = p$pairing_scheme,
    swedish_factor = p$swedish_factor,
    hunt_factor = p$hunt_factor,
    hunt_raw_sign_multiplier = p$hunt_raw_sign_multiplier,
    trait_order = load_chr$trait_order,
    species = load_chr$species,
    swedish_abs_rank = sw_rank,
    hunt_abs_rank = hu_rank,
    swedish_loading_raw = sw_load,
    hunt_loading_raw = hu_load,
    hunt_loading_frozen_aligned = format(num(hu_load) * mult, digits = 17, scientific = FALSE, trim = TRUE),
    swedish_abs_loading = format(abs(num(sw_load)), digits = 17, scientific = FALSE, trim = TRUE),
    hunt_abs_loading = format(abs(num(hu_load)), digits = 17, scientific = FALSE, trim = TRUE),
    signed_direction_concordant_after_frozen_alignment = bool_chr(sign(num(sw_load)) == sign(num(hu_load) * mult)),
    shared_top20 = bool_chr(num(sw_rank) <= 20 & num(hu_rank) <= 20),
    shared_top40 = bool_chr(num(sw_rank) <= 40 & num(hu_rank) <= 40),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
stopifnot(nrow(rank_rows) == 964L)

# Use the row-local factor to avoid any dependency on pair order.
f2d_alignment_expected <- vapply(seq_len(nrow(rank_rows)), function(i) {
  trait_i <- as.integer(rank_rows$trait_order[i])
  if (rank_rows$hunt_factor[i] == "F1") load_num$HUNT_loading_aligned_to_Swedish_F1[trait_i] else load_num$HUNT_loading_aligned_to_Swedish_F2[trait_i]
}, numeric(1))
f2d_alignment_max_abs_diff <- max(abs(num(rank_rows$hunt_loading_frozen_aligned) - f2d_alignment_expected))
stopifnot(f2d_alignment_max_abs_diff < 1e-15)

# The frozen core-species table is content-checked against matched Top20 rows and directions.
core_primary <- core_num %>%
  filter(record_type == "shared_primary_top20_species", pair_id %in% c("F1_F1", "F2_F2")) %>%
  arrange(pair_id, species)
rank_primary <- rank_rows %>%
  filter(pair_id %in% c("F1_F1", "F2_F2"), shared_top20 == "TRUE") %>%
  arrange(pair_id, species)
stopifnot(nrow(core_primary) == 15L, nrow(rank_primary) == 15L,
          identical(paste(core_primary$pair_id, core_primary$species), paste(rank_primary$pair_id, rank_primary$species)),
          max(abs(core_primary$swedish_loading_signed - num(rank_primary$swedish_loading_raw))) < 1e-15,
          max(abs(core_primary$hunt_loading_aligned_signed - num(rank_primary$hunt_loading_frozen_aligned))) < 1e-15,
          identical(as.logical(core_primary$direction_concordant), rank_primary$signed_direction_concordant_after_frozen_alignment == "TRUE"))
core_species_records_verified <- nrow(core_primary)
write_tsv(rank_rows, file.path(out_dir, "F2D_species_rank_concordance.tsv"))

# Panel C: all 241 species for each factor pair and each requested top-k, with frozen summaries repeated.
topk_frozen <- topk_chr %>%
  filter(comparison_scope == "cross_cohort_pair",
         pair_id %in% pair_defs$pair_id,
         top_k %in% c("20", "40"))
stopifnot(nrow(topk_frozen) == 8L)

parse_common_species <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  entries <- strsplit(x, "; ", fixed = TRUE)[[1]]
  sub("\\[[0-9]+:[0-9]+\\]$", "", entries)
}

f2c_list <- lapply(seq_len(nrow(topk_frozen)), function(i) {
  z <- topk_frozen[i, ]
  base <- rank_rows %>% filter(pair_id == z$pair_id)
  k <- as.numeric(z$top_k)
  sw_in <- num(base$swedish_abs_rank) <= k
  hu_in <- num(base$hunt_abs_rank) <= k
  shared <- sw_in & hu_in
  frozen_names <- sort(parse_common_species(z$common_species_rank_pairs))
  extracted_names <- sort(base$species[shared])
  if (!identical(frozen_names, extracted_names)) {
    stop("Top-k membership mismatch for ", z$pair_id, " K=", k)
  }
  if (sum(shared) != as.numeric(z$overlap_count)) {
    stop("Top-k overlap-count mismatch for ", z$pair_id, " K=", k)
  }
  base %>%
    mutate(
      top_k = z$top_k,
      swedish_in_topk = bool_chr(sw_in),
      hunt_in_topk = bool_chr(hu_in),
      shared_species = bool_chr(shared),
      in_union = bool_chr(sw_in | hu_in),
      frozen_overlap_count = z$overlap_count,
      extracted_overlap_count = as.character(sum(shared)),
      frozen_jaccard = z$jaccard,
      frozen_overlap_fraction = z$overlap_fraction,
      frozen_uniform_random_expected_overlap = z$uniform_random_expected_overlap,
      frozen_uniform_random_enrichment_descriptive = z$uniform_random_enrichment_descriptive,
      frozen_common_direction_concordance = z$common_direction_concordance,
      frozen_pointwise_empirical_p_descriptive = z$pointwise_empirical_p_descriptive,
      frozen_propensity_sensitivity_pointwise_empirical_p = z$propensity_sensitivity_pointwise_empirical_p,
      frozen_common_species_rank_pairs = z$common_species_rank_pairs,
      membership_verification = "EXACT_MATCH_TO_FROZEN_TOPK_ROW"
    )
})
f2c <- bind_rows(f2c_list) %>%
  arrange(match(pair_id, pair_defs$pair_id), num(top_k), num(trait_order))
stopifnot(nrow(f2c) == 1928L)
write_tsv(f2c, file.path(out_dir, "F2C_topk_overlap.tsv"))

# Final workbook S6A/S6B are checked record-by-record against the frozen Phase7 TSVs.
s6a <- as.data.frame(wb[["S6A_RBO_Permutation"]], stringsAsFactors = FALSE)
names(s6a) <- names(rbo_num)
s6a <- s6a %>% filter(comparison_scope %in% unique(rbo_num$comparison_scope))
s6a_full_check <- compare_frames(s6a, rbo_num,
                                 keys = c("comparison_scope", "pair_id", "rbo_p"), tolerance = 1e-10)
s6a_cross12 <- s6a %>% filter(comparison_scope == "cross_cohort_pair")
rbo_cross12 <- rbo_num %>% filter(comparison_scope == "cross_cohort_pair")
s6a_cross12_check <- compare_frames(s6a_cross12, rbo_cross12,
                                    keys = c("comparison_scope", "pair_id", "rbo_p"), tolerance = 1e-10)
stopifnot(s6a_cross12_check$n == 12L)

s6b <- as.data.frame(wb[["S6B_TopK_Curves"]], stringsAsFactors = FALSE)
names(s6b) <- names(topk_num)
s6b <- s6b %>% filter(comparison_scope %in% unique(topk_num$comparison_scope))
s6b_full_check <- compare_frames(s6b, topk_num,
                                 keys = c("comparison_scope", "pair_id", "top_k"), tolerance = 1e-10)
s6b_top8 <- s6b %>%
  filter(comparison_scope == "cross_cohort_pair", pair_id %in% pair_defs$pair_id,
         as.numeric(top_k) %in% c(20, 40))
topk8_workbook_check <- compare_frames(s6b_top8, topk_num %>%
                                         filter(comparison_scope == "cross_cohort_pair",
                                                pair_id %in% pair_defs$pair_id,
                                                top_k %in% c(20, 40)),
                                       keys = c("comparison_scope", "pair_id", "top_k"), tolerance = 1e-10)
stopifnot(topk8_workbook_check$n == 8L)

# Final S5A provides the matched-pair Top20/Top40 summary and direction concordance.
s5a <- as.data.frame(wb[["S5A_TopK_Overlap"]], stringsAsFactors = FALSE)
names(s5a) <- c("factor", "top_rule", "source_top_n", "observed_overlap",
                "random_overlap_expected", "overlap_enrichment",
                "overlap_direction_concordance", "exact_hypergeometric_P")
s5a <- s5a %>%
  mutate(factor_numeric = suppressWarnings(as.numeric(factor))) %>%
  filter(is.finite(factor_numeric), source_top_n %in% c(20, 40)) %>%
  mutate(pair_id = paste0("F", factor_numeric, "_F", factor_numeric), top_k = as.numeric(source_top_n)) %>%
  arrange(pair_id, top_k)
matched_topk_ref <- topk_num %>%
  filter(comparison_scope == "cross_cohort_pair", pair_id %in% c("F1_F1", "F2_F2"),
         top_k %in% c(20, 40)) %>% arrange(pair_id, top_k)
stopifnot(nrow(s5a) == 4L,
          identical(paste(s5a$pair_id, s5a$top_k), paste(matched_topk_ref$pair_id, matched_topk_ref$top_k)),
          all(as.numeric(s5a$observed_overlap) == as.numeric(matched_topk_ref$overlap_count)),
          max(abs(as.numeric(s5a$random_overlap_expected) - matched_topk_ref$uniform_random_expected_overlap)) < 1e-9,
          max(abs(as.numeric(s5a$overlap_enrichment) - matched_topk_ref$uniform_random_enrichment_descriptive)) < 1e-9,
          max(abs(as.numeric(s5a$overlap_direction_concordance) - matched_topk_ref$common_direction_concordance)) < 1e-12)
matched_direction_concordance_max_abs_diff <- max(abs(as.numeric(s5a$overlap_direction_concordance) - matched_topk_ref$common_direction_concordance))

# Final S5B is checked at the shared-species level, including ranks, aligned loadings and direction.
s5b <- as.data.frame(wb[["S5B_Shared_Species"]], stringsAsFactors = FALSE)
names(s5b) <- c("factor_pair", "top_size", "trait_order", "species",
                "swedish_accession", "hunt_accession", "swedish_loading",
                "aligned_hunt_loading", "same_direction", "swedish_abs_rank", "hunt_abs_rank")
s5b <- s5b %>%
  mutate(factor_numeric = suppressWarnings(as.numeric(factor_pair))) %>%
  filter(is.finite(factor_numeric), top_size %in% c(20, 40)) %>%
  mutate(pair_id = paste0("F", factor_numeric, "_F", factor_numeric),
         top_k = as.numeric(top_size)) %>%
  arrange(pair_id, top_k, species)
s5b_expected <- bind_rows(lapply(c(20, 40), function(k) {
  rank_rows %>%
    filter(pair_id %in% c("F1_F1", "F2_F2"), num(swedish_abs_rank) <= k, num(hunt_abs_rank) <= k) %>%
    mutate(top_k = k)
})) %>%
  mutate(
    swedish_accession = load_chr$swedish_accession[match(trait_order, load_chr$trait_order)],
    hunt_accession = load_chr$hunt_accession[match(trait_order, load_chr$trait_order)]
  ) %>% arrange(pair_id, top_k, species)
stopifnot(nrow(s5b) == 47L, nrow(s5b_expected) == 47L,
          identical(paste(s5b$pair_id, s5b$top_k, s5b$species),
                    paste(s5b_expected$pair_id, s5b_expected$top_k, s5b_expected$species)),
          identical(as.character(s5b$swedish_accession), s5b_expected$swedish_accession),
          identical(as.character(s5b$hunt_accession), s5b_expected$hunt_accession),
          max(abs(as.numeric(s5b$swedish_loading) - num(s5b_expected$swedish_loading_raw))) < 1e-10,
          max(abs(as.numeric(s5b$aligned_hunt_loading) - num(s5b_expected$hunt_loading_frozen_aligned))) < 1e-10,
          identical(as.logical(s5b$same_direction), s5b_expected$signed_direction_concordant_after_frozen_alignment == "TRUE"),
          all(as.numeric(s5b$swedish_abs_rank) == num(s5b_expected$swedish_abs_rank)),
          all(as.numeric(s5b$hunt_abs_rank) == num(s5b_expected$hunt_abs_rank)))
s5b_shared_species_records_verified <- nrow(s5b)

# Panel E: cross-cohort, within-cohort and global two-dimensional correspondence.
rbo_primary <- rbo_chr %>%
  filter(rbo_p == "0.95",
         (comparison_scope == "cross_cohort_pair" & pair_id %in% pair_defs$pair_id) |
         comparison_scope == "within_cohort_factor_overlap")
stopifnot(nrow(rbo_primary) == 6L)

f2e_pairs <- rbo_primary %>%
  transmute(
    record_type = ifelse(comparison_scope == "cross_cohort_pair", "factor_pair_full_profile", "within_cohort_factor_pair"),
    comparison_scope, pair_id, pair_status, pairing_scheme,
    cohort_A = ifelse(comparison_scope == "cross_cohort_pair", "Swedish", sub("_F1_F2$", "", pair_id)),
    factor_A = ifelse(comparison_scope == "cross_cohort_pair", swedish_factor, "F1"),
    cohort_B = ifelse(comparison_scope == "cross_cohort_pair", "HUNT", sub("_F1_F2$", "", pair_id)),
    factor_B = ifelse(comparison_scope == "cross_cohort_pair", hunt_factor, "F2"),
    hunt_raw_sign_multiplier,
    signed_loading_definition = ifelse(
      comparison_scope == "cross_cohort_pair",
      "Swedish raw loading versus HUNT raw loading multiplied by frozen HUNT factor sign multiplier",
      "Raw F1 versus raw F2 within the named cohort"
    ),
    n_species = list_length,
    pearson_signed = pearson_signed_all241,
    spearman_signed = spearman_signed_all241,
    pearson_absolute = pearson_absolute_all241,
    spearman_absolute = spearman_absolute_all241,
    direction_concordance = direction_concordance_all241,
    one_dim_principal_angle_deg,
    global_2d_principal_angle_1_deg = "",
    global_2d_principal_angle_2_deg = "",
    exact_comparison_definition = ifelse(
      comparison_scope == "cross_cohort_pair",
      paste0("All 241 species: Swedish ", swedish_factor, " versus sign-oriented HUNT ", hunt_factor),
      paste0("All 241 species: ", sub("_F1_F2$", "", pair_id), " raw F1 versus raw F2")
    ),
    evidence_role = ifelse(pair_status == "primary_frozen", "MATCHED_FACTOR", ifelse(pair_status == "off_diagonal_sensitivity", "CROSS_FACTOR", "WITHIN_COHORT_DIAGNOSTIC")),
    source_field = "GLOBAL_RBO_ALL_PAIRS.tsv; p=0.95 row"
  )

angle_parts <- strsplit(rbo_primary$global_2d_principal_angles_deg[which(rbo_primary$comparison_scope == "cross_cohort_pair")[1]], ";", fixed = TRUE)[[1]]
stopifnot(length(angle_parts) == 2L, identical(angle_parts, c("58.955", "80.892")))
f2e_subspace <- data.frame(
  record_type = "global_2d_loading_subspace",
  comparison_scope = "cross_cohort_2d_subspace",
  pair_id = "Swedish_F1_F2_vs_HUNT_F1_F2",
  pair_status = "global_diagnostic",
  pairing_scheme = "identity_frozen",
  cohort_A = "Swedish", factor_A = "F1+F2 subspace",
  cohort_B = "HUNT", factor_B = "F1+F2 subspace",
  hunt_raw_sign_multiplier = "not_applicable_rotation_invariant",
  signed_loading_definition = "Principal angles between the complete two-dimensional loading subspaces; invariant to factor signs",
  n_species = "241",
  pearson_signed = "", spearman_signed = "", pearson_absolute = "", spearman_absolute = "",
  direction_concordance = "", one_dim_principal_angle_deg = "",
  global_2d_principal_angle_1_deg = angle_parts[1],
  global_2d_principal_angle_2_deg = angle_parts[2],
  exact_comparison_definition = "Column space of Swedish 241x2 loading matrix versus column space of HUNT 241x2 loading matrix",
  evidence_role = "SUBSPACE_DIAGNOSTIC",
  source_field = "GLOBAL_RBO_ALL_PAIRS.tsv: global_2d_principal_angles_deg",
  stringsAsFactors = FALSE
)
f2e <- bind_rows(f2e_pairs, f2e_subspace)
stopifnot(nrow(f2e) == 7L)
write_tsv(f2e, file.path(out_dir, "F2E_factor_correspondence.tsv"))

# Cross-check exact raw correlations against final S4A and matched alignment against S4B.
s4a <- wb[["S4A_Factor_Correlations"]]
names(s4a) <- c("layer", "metric", "swedish_factor", "hunt_factor", "value")
s4a_species <- s4a %>%
  filter(layer == "species_loadings", metric %in% c("pearson", "spearman")) %>%
  mutate(swedish_factor = paste0("F", swedish_factor), hunt_factor = paste0("F", hunt_factor), value = as.numeric(value))
cross95 <- rbo_num %>%
  filter(comparison_scope == "cross_cohort_pair", abs(rbo_p - 0.95) < 1e-12) %>%
  select(swedish_factor, hunt_factor, hunt_raw_sign_multiplier, pearson_signed_all241, spearman_signed_all241)
cor_check <- cross95 %>%
  left_join(s4a_species %>% filter(metric == "pearson") %>% select(swedish_factor, hunt_factor, raw_pearson = value), by = c("swedish_factor", "hunt_factor")) %>%
  left_join(s4a_species %>% filter(metric == "spearman") %>% select(swedish_factor, hunt_factor, raw_spearman = value), by = c("swedish_factor", "hunt_factor")) %>%
  mutate(pearson_diff = pearson_signed_all241 - raw_pearson * hunt_raw_sign_multiplier,
         spearman_diff = spearman_signed_all241 - raw_spearman * hunt_raw_sign_multiplier)
correlation_crosscheck_max_abs_diff <- max(abs(c(cor_check$pearson_diff, cor_check$spearman_diff)))
stopifnot(correlation_crosscheck_max_abs_diff < 1e-9)

s4b_check <- s4b %>%
  left_join(cross95, by = c("swedish_factor", "hunt_factor")) %>%
  mutate(
    sign_diff = sign_applied_to_hunt - hunt_raw_sign_multiplier,
    pearson_aligned_diff = as.numeric(species_pearson_aligned) - pearson_signed_all241,
    spearman_aligned_diff = as.numeric(species_spearman_aligned) - spearman_signed_all241
  )
s4b_alignment_max_abs_diff <- max(abs(c(
  s4b_check$sign_diff, s4b_check$pearson_aligned_diff,
  s4b_check$spearman_aligned_diff
)))
stopifnot(nrow(s4b_check) == 2L, s4b_alignment_max_abs_diff < 1e-9)

# Verify all six pair/within-cohort correspondence rows against S6A at p=0.95,
# plus the seventh global subspace row against the frozen principal-angle field.
s6a_primary <- s6a %>%
  filter(abs(as.numeric(rbo_p) - 0.95) < 1e-12,
         (comparison_scope == "cross_cohort_pair" & pair_id %in% pair_defs$pair_id) |
           comparison_scope == "within_cohort_factor_overlap") %>%
  arrange(match(pair_id, f2e_pairs$pair_id))
f2e_pairs_check <- f2e_pairs %>% arrange(match(pair_id, s6a_primary$pair_id))
stopifnot(nrow(s6a_primary) == 6L, nrow(f2e_pairs_check) == 6L,
          identical(as.character(s6a_primary$pair_id), as.character(f2e_pairs_check$pair_id)))
factor_correspondence_max_abs_diff <- max(abs(c(
  as.numeric(s6a_primary$pearson_signed_all241) - as.numeric(f2e_pairs_check$pearson_signed),
  as.numeric(s6a_primary$spearman_signed_all241) - as.numeric(f2e_pairs_check$spearman_signed),
  as.numeric(s6a_primary$pearson_absolute_all241) - as.numeric(f2e_pairs_check$pearson_absolute),
  as.numeric(s6a_primary$spearman_absolute_all241) - as.numeric(f2e_pairs_check$spearman_absolute),
  as.numeric(s6a_primary$direction_concordance_all241) - as.numeric(f2e_pairs_check$direction_concordance),
  as.numeric(s6a_primary$one_dim_principal_angle_deg) - as.numeric(f2e_pairs_check$one_dim_principal_angle_deg)
)), na.rm = TRUE)
stopifnot(factor_correspondence_max_abs_diff < 1e-12,
          identical(c(f2e_subspace$global_2d_principal_angle_1_deg,
                      f2e_subspace$global_2d_principal_angle_2_deg), angle_parts))
factor_correspondence_records_verified <- nrow(f2e)

# Hard gate: this documentation correction must not alter any F2A-F2E bytes.
expected_f2_hashes <- c(
  F2A_factor_loadings.tsv = "886afb533c5b55d793e217f7796d79b850f9a8285801e19cb611ec1ee862e865",
  F2B_RBO_matrix.tsv = "a7a2147fc8436427f6791c59cee9536cddf26f02689bf65962f22987a2d975a9",
  F2C_topk_overlap.tsv = "ca08e3f6e7f2ff65e443566525b5d337a4e7d72b7b916b0859d875297cc3149e",
  F2D_species_rank_concordance.tsv = "2d8093b6a6a65e4d1c3daf6a522fea45b8961df2940eb519e9d9c7b4d6c113c4",
  F2E_factor_correspondence.tsv = "752fa9d0772ec997e1c9f315c317d39abf05c4dd64fd69c8e18bc68e91cfe052"
)
actual_f2_hashes <- vapply(file.path(out_dir, names(expected_f2_hashes)), function(x) {
  digest(file = x, algo = "sha256", serialize = FALSE)
}, character(1))
stopifnot(identical(unname(actual_f2_hashes), unname(expected_f2_hashes)))
f2_outputs_byte_identical <- TRUE

# Input-only SHA-256 manifest. Each role states what the script actually reads.
input_roles <- c(
  "CONTENT_CHECKED: primary frozen 241-species loading/rank table; source of F2A/F2C/F2D",
  "CONTENT_CHECKED: all frozen Phase7 RBO rows; source of F2B/F2E and S6A comparison",
  "CONTENT_CHECKED: all frozen Phase7 Top-k rows and exact common-species strings; source of F2C and S6B comparison",
  "CONTENT_CHECKED: 15 matched-pair shared Top20 records, loadings and directions compared with F2D",
  "CONTENT_CHECKED: saved seeds, permutation count and p/K grids located in executed Phase7 script text",
  "CONTENT_CHECKED: primary and sensitivity seed statements located in null-calibration text",
  "CONTENT_CHECKED: SHA-256 entries for the three consumed Phase7 numerical TSVs compared with current files",
  "CONTENT_CHECKED: 241-species taxonomy identifiers and signed loading fields compared with frozen loading table",
  "CONTENT_CHECKED: seven named worksheets read; all requested numerical/record comparisons executed",
  "CONTENT_CHECKED: seven worksheet rows and frozen controlling status verified in source crosswalk",
  "CONTENT_CHECKED: four rounded Figure 2 claims located in frozen V3 manuscript text",
  "ARCHIVE_IDENTITY_ONLY: V3 heading and revision marker checked; no Figure 2 numerical content check claimed",
  "CONTENT_CHECKED: frozen ledger records E004-E008 compared with expected rounded values"
)
manifest <- data.frame(
  input_id = names(src),
  absolute_path = unname(src),
  role = input_roles,
  bytes = file.info(src)$size,
  modified = format(file.info(src)$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  sha256 = vapply(src, function(x) digest(file = x, algo = "sha256", serialize = FALSE), character(1)),
  stringsAsFactors = FALSE
)
write_tsv(manifest, file.path(out_dir, "INPUT_SHA256_MANIFEST.tsv"))

dict <- c(
  "# Figure 2 Redesign Data Dictionary",
  "",
  "This package is a frozen-data extraction and verification package. It does not rerun GFA, RBO, permutations, rank tests or significance calculations.",
  "",
  "## Provenance decision",
  "",
  "The project directory named `PHASE16_SCZ_THREE_WAY_SIGNAL_AUDIT` contains SCZ chr3 locus/colocalization analyses and no global GFA RBO results. The final supplementary-table crosswalk, V3 manuscript and frozen evidence ledger identify `PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/GLOBAL_RBO_ALL_PAIRS.tsv` as the final global RBO source. The Phase7 source is therefore used here; this source-name mismatch is explicitly flagged rather than silently reinterpreted.",
  "",
  "## F2A_factor_loadings.tsv",
  "",
  "One row per matched species (241 rows). `*_loading_raw` and `*_abs_rank` are copied without rounding from `BIDIRECTIONAL_SPECIES_LOADING_H2.tsv`. HUNT raw signs are never overwritten. The HUNT sign multipliers are read from final worksheet `S4B_Factor_Matching` (F1=-1, F2=+1), then used to verify—not replace—the frozen aligned columns in F2A and F2D. `*_direction_raw` is derived only from the sign of the raw loading. Blank taxonomy means no unique GTDB R232 tree endpoint was resolved, not a missing microbial phenotype.",
  "",
  "## F2B_RBO_matrix.tsv",
  "",
  "Twelve frozen rows: four Swedish-HUNT factor pairs times the saved p grid (0.90, 0.95, 0.975). `primary_matrix_p=TRUE` marks the four p=0.95 rows used for the main 2x2 matrix. All original Phase7 fields are retained, including raw information-stratified empirical P, fixed-pair and all-four-pair max-statistic P fields, the propensity-sensitivity P, ranking rule and RBO formula. Permutation count, seeds, grid and null definitions are copied from the executed Phase7 script/report.",
  "",
  "## F2C_topk_overlap.tsv",
  "",
  "All 241 species for each of the four Swedish-HUNT factor pairs at Top20 and Top40 (1,928 rows). Membership flags permit exact reconstruction of both sets, their union and intersection. `frozen_*` summary fields are copied from the corresponding Phase7 Top-k row. Membership was extracted from frozen absolute-loading ranks and required to match both the saved overlap count and the saved species string exactly; it was not inferred from a plotted image.",
  "",
  "## F2D_species_rank_concordance.tsv",
  "",
  "All 241 paired ranks for each of the four Swedish-HUNT factor comparisons (964 rows). Raw signed loadings are preserved. `hunt_loading_frozen_aligned` is an explicitly separate visualization/correspondence field produced with the saved factor-specific sign multiplier. `shared_top20` and `shared_top40` use the two frozen absolute ranks for that comparison.",
  "",
  "## F2E_factor_correspondence.tsv",
  "",
  "Seven records: four cross-cohort factor pairs, two within-cohort F1-F2 diagnostics, and one global two-dimensional subspace record. Pairwise correlations and one-dimensional angles are copied from the p=0.95 Phase7 rows. The global record contains the exact principal angles 58.955 and 80.892 degrees. Signed cross-cohort correlations compare Swedish raw loadings with HUNT loadings after the separately recorded frozen sign multiplier; absolute-loading correlations are sign invariant.",
  "",
  "## Missing-value conventions",
  "",
  "Empty fields represent either (i) unavailable GTDB mapping for an otherwise valid species, or (ii) structurally inapplicable statistics such as pairwise correlations on the global subspace record or fixed-pair family P values for off-diagonal comparisons. No numerical value was imputed.",
  "",
  "## Final-workbook verification coverage",
  "",
  "The script explicitly reads `S3_Loadings_PVE`, `S4A_Factor_Correlations`, `S4B_Factor_Matching`, `S5A_TopK_Overlap`, `S5B_Shared_Species`, `S6A_RBO_Permutation`, and `S6B_TopK_Curves`. It compares 241 species and four loading profiles; all 12 cross-cohort RBO records; all eight Top20/Top40 cross-cohort set definitions; matched-pair direction concordance and shared-species records; and all seven factor-correspondence records. Workbook display rounding is tolerated only within documented numerical tolerances; exported F2A-F2E bytes are protected by fixed SHA-256 gates.",
  ""
)
writeLines(dict, file.path(out_dir, "FIGURE2_DATA_DICTIONARY.md"), useBytes = TRUE)

missing_cells <- function(df) sum(is.na(df) | (is.character(as.matrix(df)) & as.matrix(df) == ""))
missing_cols <- function(df) {
  z <- vapply(df, function(x) sum(is.na(x) | (is.character(x) & x == "")), integer(1))
  z <- z[z > 0]
  if (!length(z)) return("none")
  paste(paste0(names(z), "=", z), collapse = "; ")
}

f2a_missing_tax <- sum(is.na(f2a$gtdb_tree_tip) | f2a$gtdb_tree_tip == "")
f2a_missing_loading <- sum(is.na(f2a[, grep("loading_raw$|abs_rank$", names(f2a))]))

qc <- c(
  "# Figure 2 Numerical QC Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Checks executed by the script",
  "",
  paste0("- PASS: all seven required final supplementary worksheets were read: ", paste(required_worksheets, collapse = ", "), "."),
  paste0("- PASS: S3 comparison covered 241 species x four loading profiles; maximum display-rounding difference = ", format(workbook_loading_max_abs_diff, scientific = TRUE, digits = 3), "; maximum rank difference = ", workbook_rank_max_abs_diff, "."),
  paste0("- PASS: HUNT sign multipliers were read from S4B (F1=", hunt_sign_map[["F1"]], ", F2=", hunt_sign_map[["F2"]], "); F2A/F2D aligned-column maximum difference = ", format(max(f2a_alignment_max_abs_diff, f2d_alignment_max_abs_diff), scientific = TRUE, digits = 3), "."),
  paste0("- PASS: S6A full-table records = ", s6a_full_check$n, "; requested cross-cohort RBO records = ", s6a_cross12_check$n, "; maximum numerical difference = ", format(max(s6a_full_check$max_numeric_diff, s6a_cross12_check$max_numeric_diff), scientific = TRUE, digits = 3), "."),
  paste0("- PASS: S6B full-table records = ", s6b_full_check$n, "; requested Top20/Top40 cross-cohort records = ", topk8_workbook_check$n, "; maximum numerical difference = ", format(max(s6b_full_check$max_numeric_diff, topk8_workbook_check$max_numeric_diff), scientific = TRUE, digits = 3), "."),
  paste0("- PASS: all eight Top20/Top40 cross-cohort membership sets exactly matched saved species strings and counts: ", paste(paste0(topk_frozen$pair_id, "@", topk_frozen$top_k, "=", topk_frozen$overlap_count), collapse = "; "), "."),
  paste0("- PASS: S5A matched-factor direction concordance maximum difference = ", format(matched_direction_concordance_max_abs_diff, scientific = TRUE, digits = 3), "; S5B shared-species records verified = ", s5b_shared_species_records_verified, "."),
  paste0("- PASS: S4A aligned correlation maximum difference = ", format(correlation_crosscheck_max_abs_diff, scientific = TRUE, digits = 3), "; S4B aligned-field maximum difference = ", format(s4b_alignment_max_abs_diff, scientific = TRUE, digits = 3), "."),
  paste0("- PASS: all ", factor_correspondence_records_verified, " F2E factor-correspondence records were verified; maximum pairwise-field difference = ", format(factor_correspondence_max_abs_diff, scientific = TRUE, digits = 3), "."),
  paste0("- PASS: frozen core-species content records verified = ", core_species_records_verified, "; Phase7 hash-ledger entries verified = ", length(phase7_hash_targets), "."),
  "- PASS: V3 manuscript checks were limited to four rounded claims; evidence-ledger checks were limited to records E004-E008; V3 handoff was checked only for archive identity.",
  paste0("- PASS: F2A rows=", nrow(f2a), ", F2B rows=", nrow(f2b), ", F2C rows=", nrow(f2c), ", F2D rows=", nrow(rank_rows), ", F2E rows=", nrow(f2e), "."),
  paste0("- PASS: F2A unresolved GTDB tree tips=", f2a_missing_tax, "; missing loading/rank values=", f2a_missing_loading, ". Structural blanks in other files: F2B {", missing_cols(f2b), "}; F2C {", missing_cols(f2c), "}; F2D {", missing_cols(rank_rows), "}; F2E {", missing_cols(f2e), "}."),
  "- PASS: fixed pre-correction SHA-256 gates confirm F2A-F2E are byte-identical after this documentation correction.",
  "",
  "## Independent manual cross-checks",
  "",
  "- The seven worksheet names, their displayed headers and their expected row families were visually inspected during this correction.",
  "- The generated TSV headers and report section labels were manually inspected for naming and documentation consistency.",
  "- These manual checks are not the basis for the numerical PASS statements above; those statements are generated from executable checks in this script.",
  "",
  "## Checks not performed",
  "",
  "- No GFA model was fitted or rerun.",
  "- No RBO statistic, permutation distribution, rank, P value, adjusted P value or significance result was recomputed.",
  "- No loading sign was changed and no species identifier was renamed.",
  "- No content-wide audit of the V3 handoff was performed; only its archive identity markers were checked.",
  "- No claims outside the explicitly listed workbook sheets, V3 manuscript phrases and evidence-ledger records E004-E008 were evaluated.",
  "",
  "## Source-label note",
  "",
  "The requested 'final frozen Phase16 RBO' does not exist as a global-GFA result: Phase16 is an SCZ chr3 audit. The final supplementary crosswalk and frozen manuscript identify Phase7 as the authoritative global RBO source; this mismatch is documented rather than silently relabeled.",
  ""
)
writeLines(qc, file.path(out_dir, "FIGURE2_NUMERICAL_QC_REPORT.md"), useBytes = TRUE)

# Output manifest covers the five frozen numerical exports and the three updated
# documentation/reproducibility artifacts. It intentionally excludes itself.
output_files <- c(
  "F2A_factor_loadings.tsv", "F2B_RBO_matrix.tsv", "F2C_topk_overlap.tsv",
  "F2D_species_rank_concordance.tsv", "F2E_factor_correspondence.tsv",
  "FIGURE2_DATA_DICTIONARY.md", "FIGURE2_NUMERICAL_QC_REPORT.md",
  "build_figure2_redesign_data.R"
)
output_paths <- file.path(out_dir, output_files)
stopifnot(all(file.exists(output_paths)))
output_manifest <- data.frame(
  artifact = output_files,
  bytes = as.numeric(file.info(output_paths)$size),
  sha256 = vapply(output_paths, function(x) digest(file = x, algo = "sha256", serialize = FALSE), character(1)),
  stringsAsFactors = FALSE
)
write_tsv(output_manifest, file.path(out_dir, "OUTPUT_SHA256_MANIFEST.tsv"))

message("Figure 2 redesign data package completed: ", out_dir)
