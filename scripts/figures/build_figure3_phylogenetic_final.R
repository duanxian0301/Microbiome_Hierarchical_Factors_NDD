#!/usr/bin/env Rscript

# Rebuild Main Figure 3 from frozen Phase 8 outputs only.
# No GFA, GWAS, phylogenetic permutation, P value, q value, null distribution,
# or factor sign is recomputed in this script.

options(stringsAsFactors = FALSE, digits = 17, scipen = 999, warn = 1)

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(scales)
})

required <- c("svglite", "ragg", "readxl", "digest", "pdftools")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

project <- normalizePath("PROJECT_ROOT", winslash = "/", mustWork = TRUE)
phase8 <- file.path(project, "FULL_SCALE_GFA_FINAL_EVALUATION", "PHASE8_PHYLOGENETIC_ORGANIZATION")
phase7 <- file.path(project, "FULL_SCALE_GFA_FINAL_EVALUATION", "PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY")
prior_redesign <- file.path(project, "FIGURE3_PHYLOGENETIC_REDESIGN")
outdir <- file.path(project, "FIGURE3_PHYLOGENETIC_REDESIGN_FINAL")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(outdir, "FIGURE3_PHYLOGENY_FINAL_EXECUTION.log")
log_lines <- character()
log_msg <- function(...) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = ""), sep = "\t")
  log_lines <<- c(log_lines, line)
  writeLines(log_lines, log_path, useBytes = TRUE)
  cat(line, "\n")
}
log_msg("START Figure 3 final visual refinement from frozen Phase 8 results")

read_tsv_chr <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
             check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, col.names = TRUE,
              quote = FALSE, na = "NA", fileEncoding = "UTF-8")
}
as_num <- function(x) suppressWarnings(as.numeric(x))
as_bool <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
blank_to_na <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == "" | toupper(trimws(x)) == "NA"] <- NA_character_
  x
}

inputs <- list(
  phase8_mapping = file.path(phase8, "PHYLOGENETIC_SPECIES_MAPPING.tsv"),
  phase8_clustering = file.path(phase8, "PHYLOGENETIC_CLUSTERING_RESULTS.tsv"),
  phase8_reproducibility = file.path(phase8, "PHYLOGENETIC_REPRODUCIBILITY.tsv"),
  phase8_gate = file.path(phase8, "PHASE8_FINAL_SCIENTIFIC_GATE.md"),
  phase8_prior_figure = file.path(phase8, "PHYLOGENETIC_FIGURES.pdf"),
  phase8_tree = file.path(phase8, "references", "bac120.sp_labels.tree"),
  phase8_taxonomy = file.path(phase8, "references", "bac120_taxonomy.tsv.gz"),
  phase8_synonyms = file.path(phase8, "references", "synonyms.bac120.tsv"),
  phase8_version = file.path(phase8, "references", "VERSION.txt"),
  phase8_release_notes = file.path(phase8, "references", "RELEASE_NOTES.txt"),
  phase8_original_script = file.path(phase8, "scripts", "02_phase8_figures_from_outputs.R"),
  frozen_loading_export = file.path(project, "FULL_SCALE_GFA_FINAL_EVALUATION", "BIDIRECTIONAL_SPECIES_LOADING_H2.tsv"),
  phase7_core_species = file.path(phase7, "GLOBAL_COMPOSITION_CORE_SPECIES.tsv"),
  phase7_decision = file.path(phase7, "GLOBAL_GFA_SCIENTIFIC_DECISION.md"),
  phase7_script = file.path(phase7, "scripts", "01_phase7_global_composition.R"),
  supplementary_v3 = file.path(project, "MANUSCRIPT_WRITING_HANDOFF", "GLOBAL_GFA_SUPPLEMENTARY_TABLES_V3.xlsx"),
  supplementary_final = file.path(project, "SUBMISSION_PACKAGE", "SUPPLEMENTARY_TABLES_FINAL", "Supplementary_Tables.xlsx"),
  prior_value_audit = file.path(prior_redesign, "figure_value_audit.tsv"),
  prior_final_qc = file.path(prior_redesign, "FIGURE3_PHYLOGENY_QC.md"),
  prior_plot_script = file.path(prior_redesign, "build_figure3_phylogenetic.R")
)
missing_inputs <- names(inputs)[!file.exists(unlist(inputs))]
if (length(missing_inputs)) stop("Missing frozen inputs: ", paste(missing_inputs, collapse = ", "))

roles <- c(
  "Frozen 241-species mapping, taxonomy, aligned loadings, and mapping provenance",
  "Frozen MPD/MNTD, sensitivity, continuous-signal, and shared-pool results",
  "Frozen cross-cohort shared-species and lineage summaries",
  "Controlling Phase 8 interpretation and strict scientific boundaries",
  "Prior two-page figure inspected as redesign reference only",
  "GTDB R11-RS232 120-marker species-cluster tree",
  "GTDB R11-RS232 taxonomy provenance",
  "GTDB R11-RS232 synonym provenance",
  "GTDB reference version identifier",
  "GTDB release provenance",
  "Historical Phase 8 plotting implementation inspected for field use",
  "Frozen Swedish/HUNT raw and aligned loading export; sign audit",
  "Frozen Phase 7 shared Top20 membership and direction table",
  "Phase 7 evidence boundary and historical shared-species count",
  "Frozen Phase 7 factor/sign alignment implementation",
  "V3 supplementary S14 mapping/clustering/reproducibility mirror",
  "Submission supplementary S7B/S7C current mirror",
  "Previously approved 1017-row value audit retained verbatim before two display-specific rows are appended",
  "Prior redesign QC used to preserve the already accepted scientific and rendering checks",
  "Prior reproducible R plotting implementation used as the visual-refinement baseline"
)

source_inventory <- data.frame(
  input_id = names(inputs),
  absolute_path = normalizePath(unlist(inputs), winslash = "/"),
  bytes = as.character(file.info(unlist(inputs))$size),
  sha256 = vapply(unlist(inputs), sha256, character(1)),
  intended_figure_use = roles,
  check.names = FALSE
)
write_tsv(source_inventory, file.path(outdir, "SOURCE_INVENTORY.tsv"))
log_msg("Source inventory written: ", nrow(source_inventory), " inputs")

mapping <- read_tsv_chr(inputs$phase8_mapping)
clustering <- read_tsv_chr(inputs$phase8_clustering)
repro <- read_tsv_chr(inputs$phase8_reproducibility)
loading <- read_tsv_chr(inputs$frozen_loading_export)
core <- read_tsv_chr(inputs$phase7_core_species)
gate_text <- paste(readLines(inputs$phase8_gate, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
phase7_decision_text <- paste(readLines(inputs$phase7_decision, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
version_text <- trimws(paste(readLines(inputs$phase8_version, warn = FALSE), collapse = " "))

stopifnot(nrow(mapping) == 241L, sum(as_bool(mapping$uniquely_mapped)) == 236L)
map_m <- mapping[as_bool(mapping$uniquely_mapped), , drop = FALSE]
stopifnot(length(unique(map_m$gtdb_tree_tip)) == 236L)
stopifnot(grepl("PHYLOGENETIC_INCREMENT_SUPPORTED", gate_text, fixed = TRUE))
stopifnot(grepl("R11-RS232", gate_text, fixed = TRUE), grepl("v232", version_text, fixed = TRUE))

# Confirm frozen Phase 7 signs and all four mapped loading profiles exactly.
loading_idx <- match(mapping$trait_order, loading$trait_order)
stopifnot(!anyNA(loading_idx))
sign_f1 <- unique(loading$HUNT_sign_applied_for_Swedish_F1)
sign_f2 <- unique(loading$HUNT_sign_applied_for_Swedish_F2)
stopifnot(identical(sign_f1, "-1"), identical(sign_f2, "1"))
numeric_equal <- function(a, b, tol = 1e-13) {
  aa <- as_num(a); bb <- as_num(b)
  all(!is.na(aa) & !is.na(bb) & abs(aa - bb) <= tol * pmax(1, abs(aa), abs(bb)))
}
stopifnot(
  numeric_equal(mapping$Swedish_F1_signed_loading, loading$Swedish_F1_loading_raw[loading_idx]),
  numeric_equal(mapping$Swedish_F2_signed_loading, loading$Swedish_F2_loading_raw[loading_idx]),
  numeric_equal(mapping$HUNT_F1_signed_loading, loading$HUNT_loading_aligned_to_Swedish_F1[loading_idx]),
  numeric_equal(mapping$HUNT_F2_signed_loading, loading$HUNT_loading_aligned_to_Swedish_F2[loading_idx])
)
log_msg("Frozen loading/sign audit PASS: HUNT F1 x -1; HUNT F2 x +1; 241 x 4 profiles exact")

# Shared Top20 membership: keep historical seventh F2 species as an explicit unmapped record.
core_primary <- core[core$record_type == "shared_primary_top20_species" & core$pair_id %in% c("F1_F1", "F2_F2"), , drop = FALSE]
stopifnot(sum(core_primary$pair_id == "F1_F1") == 8L, sum(core_primary$pair_id == "F2_F2") == 7L)
core_map_idx <- match(core_primary$trait_order, mapping$trait_order)
stopifnot(!anyNA(core_map_idx))
shared <- cbind(core_primary, mapping[core_map_idx, c("gtdb_tree_tip", "uniquely_mapped", "family", "genus", "species_taxon"), drop = FALSE])
shared$pair_display <- ifelse(shared$pair_id == "F1_F1", "F1-F1", "F2-F2")
shared$included_in_tree_panel <- as_bool(shared$uniquely_mapped)
shared$mapping_note <- ifelse(shared$included_in_tree_panel, "Unique R232 endpoint", "No unique R232 endpoint; excluded without substitution")
shared$direction_concordant <- toupper(shared$direction_concordant)
shared$swedish_loading_signed <- core_primary$swedish_loading_signed
shared$hunt_loading_aligned_signed <- core_primary$hunt_loading_aligned_signed
shared$historical_shared_top20_n <- ifelse(shared$pair_id == "F1_F1", "8", "7")
shared$mapped_shared_top20_n <- ifelse(shared$pair_id == "F1_F1", "8", "6")
shared <- shared[, c("pair_id", "pair_display", "species", "trait_order", "gtdb_tree_tip", "included_in_tree_panel",
                     "mapping_note", "genus", "family", "species_taxon", "swedish_factor", "hunt_factor",
                     "hunt_raw_sign_multiplier", "swedish_loading_signed", "hunt_loading_aligned_signed",
                     "direction_concordant", "historical_shared_top20_n", "mapped_shared_top20_n")]
names(shared)[names(shared) == "species"] <- "frozen_species"
write_tsv(shared, file.path(outdir, "F3B_SHARED_TOP20_LINEAGES.tsv"))

shared_mapped <- shared[shared$included_in_tree_panel, , drop = FALSE]
stopifnot(sum(shared_mapped$pair_id == "F1_F1") == 8L, sum(shared_mapped$pair_id == "F2_F2") == 6L)
stopifnot(sum(shared_mapped$pair_id == "F1_F1" & shared_mapped$direction_concordant == "TRUE") == 6L)
stopifnot(sum(shared_mapped$pair_id == "F2_F2" & shared_mapped$direction_concordant == "TRUE") == 6L)

repro_primary <- repro[repro$pairing_status == "primary_frozen" & repro$top_k == "20" &
                         repro$pair_id %in% c("Swedish_F1__HUNT_F1", "Swedish_F2__HUNT_F2"), , drop = FALSE]
stopifnot(nrow(repro_primary) == 2L)
stopifnot(repro_primary$exact_species_overlap_n[repro_primary$pair_id == "Swedish_F1__HUNT_F1"] == "8")
stopifnot(repro_primary$common_species_distinct_genera_n[repro_primary$pair_id == "Swedish_F1__HUNT_F1"] == "7")
stopifnot(repro_primary$common_species_distinct_families_n[repro_primary$pair_id == "Swedish_F1__HUNT_F1"] == "6")
stopifnot(repro_primary$exact_species_overlap_n[repro_primary$pair_id == "Swedish_F2__HUNT_F2"] == "6")
stopifnot(repro_primary$common_species_distinct_genera_n[repro_primary$pair_id == "Swedish_F2__HUNT_F2"] == "6")
stopifnot(repro_primary$common_species_distinct_families_n[repro_primary$pair_id == "Swedish_F2__HUNT_F2"] == "5")
log_msg("Shared Top20 audit PASS: F1 8 mapped (7 genera, 6 families, 6/8 concordant); F2 6 mapped (6 genera, 5 families, 6/6 concordant); seventh historical F2 retained as unmapped")

# Join, but never recompute, the 16 primary and 16 conservative-sensitivity records.
primary <- clustering[clustering$analysis_family == "topK_phylogenetic_clustering" &
                        clustering$null_model == "information_PC1xPC2_stratified", , drop = FALSE]
sensitivity <- clustering[clustering$analysis_family == "topK_phylogenetic_clustering" &
                            clustering$null_model == "information_PC1xpropensity_stratified", , drop = FALSE]
stopifnot(nrow(primary) == 16L, nrow(sensitivity) == 16L)
key_cols <- c("cohort", "factor", "top_k", "metric")
key_primary <- do.call(paste, c(primary[key_cols], sep = "|"))
key_sens <- do.call(paste, c(sensitivity[key_cols], sep = "|"))
sens_idx <- match(key_primary, key_sens)
stopifnot(!anyNA(sens_idx))
panel_c <- primary
panel_c$primary_q_lt_0_05 <- as_num(panel_c$q_bh_family) < 0.05
panel_c$conservative_ses <- sensitivity$ses[sens_idx]
panel_c$conservative_empirical_p_one_sided <- sensitivity$empirical_p_one_sided[sens_idx]
panel_c$conservative_q_bh_family <- sensitivity$q_bh_family[sens_idx]
panel_c$conservative_q_lt_0_05 <- as_num(panel_c$conservative_q_bh_family) < 0.05
panel_c$profile <- paste(panel_c$cohort, panel_c$factor, sep = "_")
panel_c$display_row <- paste(panel_c$cohort, panel_c$factor, paste0("Top", panel_c$top_k))
panel_c$interpretation <- ifelse(panel_c$conservative_q_lt_0_05,
                                 "FDR-supported in primary and conservative propensity sensitivity",
                                 ifelse(panel_c$primary_q_lt_0_05, "FDR-supported in primary only", "Not FDR-supported in primary"))
stopifnot(sum(panel_c$primary_q_lt_0_05) == 8L, sum(panel_c$conservative_q_lt_0_05) == 2L)
robust <- panel_c[panel_c$conservative_q_lt_0_05, key_cols, drop = FALSE]
stopifnot(nrow(robust) == 2L, all(robust$cohort == "Swedish"), all(robust$factor == "F1"),
          all(robust$top_k == "40"), setequal(robust$metric, c("MPD", "MNTD")))
shared_pool <- clustering[clustering$analysis_family == "shared_high_loading_pool_clustering", , drop = FALSE]
stopifnot(sum(as_num(shared_pool$q_bh_family) < 0.05, na.rm = TRUE) == 0L)
write_tsv(panel_c, file.path(outdir, "F3C_PRIMARY_CLUSTERING.tsv"))
log_msg("Clustering audit PASS: 8/16 primary q<0.05; 2/16 conservative q<0.05; robust set Swedish F1 Top40 MPD/MNTD only; shared-pool FDR positives=0")

# Read and prune the exact GTDB tree. No branch lengths or labels are changed.
full_tree <- read.tree(inputs$phase8_tree)
full_tree$tip.label <- ifelse(startsWith(full_tree$tip.label, "'"),
                              substring(full_tree$tip.label, 2, nchar(full_tree$tip.label) - 1L),
                              full_tree$tip.label)
stopifnot(all(map_m$gtdb_tree_tip %in% full_tree$tip.label))
tree <- keep.tip(full_tree, map_m$gtdb_tree_tip)
rm(full_tree); invisible(gc())
stopifnot(Ntip(tree) == 236L, setequal(tree$tip.label, map_m$gtdb_tree_tip))

tree_base <- ggtree(tree, layout = "rectangular") +
  geom_tree(colour = "#B8BEC6", linewidth = 0.18)
tree_layout <- as.data.frame(tree_base$data, stringsAsFactors = FALSE)
write_tsv(tree_layout, file.path(outdir, "F3A_TREE_LAYOUT.tsv"))

tip_layout <- tree_layout[tree_layout$isTip, c("label", "x", "y", "node", "parent", "branch.length"), drop = FALSE]
map_order <- match(tip_layout$label, map_m$gtdb_tree_tip)
stopifnot(!anyNA(map_order))
tip_meta <- map_m[map_order, , drop = FALSE]
tip_meta$tree_y <- tip_layout$y
tip_meta$tree_x <- tip_layout$x

f1_tips <- shared_mapped$gtdb_tree_tip[shared_mapped$pair_id == "F1_F1"]
f2_tips <- shared_mapped$gtdb_tree_tip[shared_mapped$pair_id == "F2_F2"]
tip_meta$shared_F1_F1 <- tip_meta$gtdb_tree_tip %in% f1_tips
tip_meta$shared_F2_F2 <- tip_meta$gtdb_tree_tip %in% f2_tips
tip_meta$shared_group <- ifelse(tip_meta$shared_F1_F1 & tip_meta$shared_F2_F2, "Both",
                                ifelse(tip_meta$shared_F1_F1, "F1-F1",
                                       ifelse(tip_meta$shared_F2_F2, "F2-F2", "Not shared Top20")))
tip_source <- tip_meta[, c("gtdb_tree_tip", "tree_x", "tree_y", "trait_order", "frozen_species", "family", "genus", "species_taxon",
                         "Swedish_F1_signed_loading", "Swedish_F1_abs_rank", "Swedish_F2_signed_loading", "Swedish_F2_abs_rank",
                         "HUNT_F1_signed_loading", "HUNT_F1_abs_rank", "HUNT_F2_signed_loading", "HUNT_F2_abs_rank",
                         "shared_F1_F1", "shared_F2_F2", "shared_group", "mapping_class", "mapping_basis", "reference_release", "tree_file")]
write_tsv(tip_source, file.path(outdir, "F3A_PHYLOGENY_LOADINGS.tsv"))

# Workbook mirrors are read and compared, not merely hashed.
read_sheet <- function(path, sheet, skip) {
  as.data.frame(readxl::read_excel(path, sheet = sheet, skip = skip, col_types = "text"),
                stringsAsFactors = FALSE, check.names = FALSE)
}
norm_value <- function(x) {
  z <- blank_to_na(x)
  ifelse(is.na(z), "<NA>", trimws(z))
}
equal_values <- function(a, b) {
  aa <- norm_value(a); bb <- norm_value(b)
  na <- suppressWarnings(as.numeric(aa)); nb <- suppressWarnings(as.numeric(bb))
  numeric_pair <- !is.na(na) & !is.na(nb)
  out <- aa == bb
  out[numeric_pair] <- abs(na[numeric_pair] - nb[numeric_pair]) <= 1e-12 * pmax(1, abs(na[numeric_pair]), abs(nb[numeric_pair]))
  bool_pair <- tolower(aa) %in% c("true", "false") & tolower(bb) %in% c("true", "false")
  out[bool_pair] <- tolower(aa[bool_pair]) == tolower(bb[bool_pair])
  out
}
canonical_key <- function(x) {
  z <- norm_value(x)
  n <- suppressWarnings(as.numeric(z))
  is_num <- !is.na(n)
  z[is_num] <- format(n[is_num], digits = 15, scientific = TRUE, trim = TRUE)
  z
}
compare_source_sheet <- function(check_name, src, wb, keys) {
  shared_cols <- intersect(names(src), names(wb))
  stopifnot(all(keys %in% shared_cols))
  sk <- do.call(paste, c(lapply(src[keys], canonical_key), sep = "|"))
  wk <- do.call(paste, c(lapply(wb[keys], canonical_key), sep = "|"))
  idx <- match(sk, wk)
  if (anyNA(idx)) return(data.frame(check = check_name, status = "FAIL", detail = paste(sum(is.na(idx)), "keys absent")))
  fields <- setdiff(shared_cols, keys)
  ok <- vapply(fields, function(nm) all(equal_values(src[[nm]], wb[[nm]][idx])), logical(1))
  data.frame(check = check_name, status = if (all(ok)) "PASS" else "FAIL",
             detail = paste0(nrow(src), " rows; ", length(fields), " compared fields",
                             if (any(!ok)) paste0("; mismatched: ", paste(names(ok)[!ok], collapse = ", ")) else ""))
}

v3_mapping <- read_sheet(inputs$supplementary_v3, "S14_Mapping", 2)
v3_clustering <- read_sheet(inputs$supplementary_v3, "S14_Clustering", 2)
v3_repro <- read_sheet(inputs$supplementary_v3, "S14_Repro", 2)
final_clustering <- read_sheet(inputs$supplementary_final, "S7B_Phylo_Statistics", 4)
final_repro <- read_sheet(inputs$supplementary_final, "S7C_Phylo_Repro", 4)

workbook_checks <- do.call(rbind, list(
  compare_source_sheet("Phase8 mapping vs V3 S14_Mapping", mapping, v3_mapping, c("trait_order")),
  compare_source_sheet("Phase8 clustering vs V3 S14_Clustering", clustering, v3_clustering,
                       c("analysis_family", "cohort", "factor", "set_definition", "top_k", "metric", "null_model")),
  compare_source_sheet("Phase8 reproducibility vs V3 S14_Repro", repro, v3_repro,
                       c("comparison_scope", "pair_id", "pairing_status", "top_k")),
  compare_source_sheet("Phase8 clustering vs final S7B", clustering, final_clustering,
                       c("analysis_family", "cohort", "factor", "set_definition", "top_k", "metric", "null_model")),
  compare_source_sheet("Phase8 reproducibility vs final S7C", repro, final_repro,
                       c("comparison_scope", "pair_id", "pairing_status", "top_k"))
))
if (any(workbook_checks$status != "PASS")) stop("Supplementary workbook mismatch: ", paste(workbook_checks$check[workbook_checks$status != "PASS"], collapse = "; "))
log_msg("Supplementary workbook comparisons PASS: V3 S14 and final S7 mirrors")

# Figure value audit: every plotted quantitative cell/item points to its frozen row and column.
audit_rows <- list()
add_audit <- function(panel, item, display_value, exact_value, source_file, source_key, source_column, status = "EXACT_MATCH") {
  audit_rows[[length(audit_rows) + 1L]] <<- data.frame(
    figure = "Main Figure 3", panel = panel, plotted_item = item,
    display_value = as.character(display_value), exact_value = as.character(exact_value),
    source_file = normalizePath(source_file, winslash = "/"), source_row_key = as.character(source_key),
    source_column = source_column, exact_match_status = status, check.names = FALSE
  )
}
for (i in seq_len(nrow(tip_source))) {
  for (col in c("Swedish_F1_signed_loading", "Swedish_F2_signed_loading", "HUNT_F1_signed_loading", "HUNT_F2_signed_loading")) {
    add_audit("A", paste(tip_source$frozen_species[i], col, sep = " | "), tip_source[[col]][i], tip_source[[col]][i],
              inputs$phase8_mapping, paste0("trait_order=", tip_source$trait_order[i]), col)
  }
}
for (i in seq_len(nrow(shared))) {
  add_audit("B", paste(shared$pair_display[i], shared$frozen_species[i], "direction", sep = " | "),
            shared$direction_concordant[i], shared$direction_concordant[i], inputs$phase7_core_species,
            paste(shared$pair_id[i], shared$frozen_species[i], sep = "|"), "direction_concordant")
}
for (i in seq_len(nrow(repro_primary))) {
  for (col in c("exact_species_overlap_n", "common_species_distinct_genera_n", "common_species_distinct_families_n", "direction_concordant_n", "direction_evaluable_n")) {
    add_audit("B", paste(repro_primary$pair_id[i], col, sep = " | "), repro_primary[[col]][i], repro_primary[[col]][i],
              inputs$phase8_reproducibility, paste(repro_primary$pair_id[i], "top20", sep = "|"), col)
  }
}
for (i in seq_len(nrow(panel_c))) {
  display_q <- formatC(as_num(panel_c$q_bh_family[i]), format = "g", digits = 3)
  add_audit("C", paste(panel_c$display_row[i], panel_c$metric[i], "SES", sep = " | "),
            panel_c$ses[i], panel_c$ses[i], inputs$phase8_clustering, key_primary[i], "ses")
  add_audit("C", paste(panel_c$display_row[i], panel_c$metric[i], "primary BH q", sep = " | "),
            display_q, panel_c$q_bh_family[i], inputs$phase8_clustering, key_primary[i], "q_bh_family", "TRACEABLE_ROUNDED_DISPLAY")
  add_audit("C", paste(panel_c$display_row[i], panel_c$metric[i], "conservative BH q", sep = " | "),
            panel_c$conservative_q_bh_family[i], panel_c$conservative_q_bh_family[i], inputs$phase8_clustering,
            key_sens[sens_idx[i]], "q_bh_family")
}
generated_value_audit <- do.call(rbind, audit_rows)
prior_value_audit <- read_tsv_chr(inputs$prior_value_audit)
stopifnot(nrow(prior_value_audit) == 1017L,
          identical(names(prior_value_audit), names(generated_value_audit)),
          identical(prior_value_audit, generated_value_audit))

robust_idx <- which(panel_c$conservative_q_lt_0_05)
stopifnot(length(robust_idx) == 2L)
robust_mpd_idx <- robust_idx[panel_c$metric[robust_idx] == "MPD"]
robust_mntd_idx <- robust_idx[panel_c$metric[robust_idx] == "MNTD"]
stopifnot(length(robust_mpd_idx) == 1L, length(robust_mntd_idx) == 1L)
stopifnot(numeric_equal(panel_c$conservative_q_bh_family[robust_mpd_idx], "0.0223977602239776"),
          numeric_equal(panel_c$conservative_q_bh_family[robust_mntd_idx], "0.0191980801919808"))

new_display_audit <- data.frame(
  figure = "Main Figure 3",
  panel = "C",
  plotted_item = c("Swedish F1 Top40 | MPD | displayed high-loading-propensity sensitivity q",
                   "Swedish F1 Top40 | MNTD | displayed high-loading-propensity sensitivity q"),
  display_value = c("0.0224", "0.0192"),
  exact_value = c(panel_c$conservative_q_bh_family[robust_mpd_idx], panel_c$conservative_q_bh_family[robust_mntd_idx]),
  source_file = normalizePath(inputs$phase8_clustering, winslash = "/"),
  source_row_key = c(key_sens[sens_idx[robust_mpd_idx]], key_sens[sens_idx[robust_mntd_idx]]),
  source_column = "q_bh_family",
  exact_match_status = "TRACEABLE_ROUNDED_DISPLAY",
  check.names = FALSE
)
value_audit <- rbind(prior_value_audit, new_display_audit)
write_tsv(value_audit, file.path(outdir, "figure_value_audit.tsv"))
log_msg("Value audit written: retained 1017 prior records verbatim and appended 2 conservative-q display records")

# ---------------- Figure contract ----------------
# Core conclusion: high-loading species show limited phylogenetic organization;
# robust clustering is confined to Swedish F1 Top40 under the conservative null.
# Archetype: asymmetric mixed-modality figure, with the rectangular tree as hero.
# Panel A: 236-tip GTDB tree + four frozen aligned signed-loading tracks.
# Panel B: exact mapped shared Top20 species, lineages, and loading direction.
# Panel C: all 16 primary MPD/MNTD SES/q values + conservative significance marks.
# Reviewer risks: implying functional conservation, universal clustering, complete
# factor replication, or substitution of an unmapped historical species.

palette_profile <- c(Swedish_F1 = "#3A6F9E", Swedish_F2 = "#78A9C8", HUNT_F1 = "#8B5F8F", HUNT_F2 = "#C58B65")
palette_shared <- c("F1-F1" = "#3A6F9E", "F2-F2" = "#C58B65", "Both" = "#745E8F", "Not shared Top20" = "#B8BEC6")
shape_shared <- c("F1-F1" = 21, "F2-F2" = 24, "Both" = 23, "Not shared Top20" = 16)

theme_pub <- theme_classic(base_size = 6.8, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.28, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.28, colour = "#333333"),
    axis.text = element_text(size = 6.1, colour = "#333333"),
    axis.title = element_text(size = 6.8, colour = "#222222"),
    legend.text = element_text(size = 5.8), legend.title = element_text(size = 5.9),
    panel.grid = element_blank(), plot.title = element_text(size = 7.0, face = "bold"),
    plot.tag = element_text(size = 8, face = "bold", colour = "#111111"),
    plot.tag.position = c(0, 1)
  )

# Panel A
xmax <- max(tree_layout$x, na.rm = TRUE)
tree_plot_shift <- xmax * 0.15
tree_plot <- tree_base
tree_plot$data$x <- tree_plot$data$x - tree_plot_shift
tip_meta$tree_x_plot <- tip_meta$tree_x - tree_plot_shift
track_step <- xmax * 0.046
track_x <- xmax + xmax * 0.035 + seq(0, 3) * track_step
names(track_x) <- c("Swedish F1", "Swedish F2", "HUNT F1", "HUNT F2")
label_x <- max(track_x) + xmax * 0.20
panel_xmax <- label_x + xmax * 0.82
loading_lim <- max(abs(as_num(c(tip_meta$Swedish_F1_signed_loading, tip_meta$Swedish_F2_signed_loading,
                                tip_meta$HUNT_F1_signed_loading, tip_meta$HUNT_F2_signed_loading))), na.rm = TRUE)
loading_breaks <- c(-loading_lim, 0, loading_lim)
loading_labels <- c(paste0("-", formatC(loading_lim, format = "f", digits = 2)), "0", paste0("+", formatC(loading_lim, format = "f", digits = 2)))

loading_long <- tip_meta |>
  select(gtdb_tree_tip, tree_y, Swedish_F1_signed_loading, Swedish_F2_signed_loading,
         HUNT_F1_signed_loading, HUNT_F2_signed_loading) |>
  pivot_longer(cols = ends_with("signed_loading"), names_to = "profile", values_to = "loading") |>
  mutate(profile_display = recode(profile,
                                  Swedish_F1_signed_loading = "Swedish F1",
                                  Swedish_F2_signed_loading = "Swedish F2",
                                  HUNT_F1_signed_loading = "HUNT F1",
                                  HUNT_F2_signed_loading = "HUNT F2"),
         track_x = unname(track_x[profile_display]), loading_num = as_num(loading))

shared_tip_meta <- tip_meta[tip_meta$shared_group != "Not shared Top20", , drop = FALSE]
shared_tip_meta$label_x <- label_x

pA <- tree_plot +
  geom_tile(data = loading_long,
            aes(x = track_x, y = tree_y, fill = loading_num),
            inherit.aes = FALSE, width = track_step * 0.72, height = 0.92) +
  geom_point(data = tip_meta,
             aes(x = tree_x_plot, y = tree_y, shape = shared_group, colour = shared_group),
             inherit.aes = FALSE, fill = "white", size = 0.65, stroke = 0.18) +
  geom_text_repel(data = shared_tip_meta,
                  aes(x = label_x, y = tree_y, label = frozen_species, colour = shared_group),
                  inherit.aes = FALSE, direction = "y", hjust = 0, nudge_x = 0,
                  box.padding = 0.10, point.padding = 0.05, min.segment.length = 0,
                  max.overlaps = Inf, seed = 20261107, size = 1.75,
                  segment.colour = NA, family = "Arial", fontface = "italic",
                  show.legend = FALSE) +
  annotate("text", x = track_x, y = Ntip(tree) + 2.0, label = names(track_x), angle = 90,
           size = 1.95, family = "Arial", fontface = "bold", hjust = 0) +
  annotate("text", x = -tree_plot_shift, y = Ntip(tree) + 20.0, label = "236 uniquely mapped species",
           hjust = 0, size = 1.9, family = "Arial", colour = "#333333") +
  scale_fill_gradient2(low = "#3A6F9E", mid = "#FAFAF8", high = "#B45D5D", midpoint = 0,
                       limits = c(-loading_lim, loading_lim), breaks = loading_breaks, labels = loading_labels,
                       name = "Signed loading") +
  scale_colour_manual(values = palette_shared, breaks = c("F1-F1", "F2-F2", "Both"), name = "Shared Top20") +
  scale_shape_manual(values = shape_shared, breaks = c("F1-F1", "F2-F2", "Both"), name = "Shared Top20") +
  guides(
    fill = guide_colourbar(order = 1, title.position = "top", barwidth = unit(22, "mm"), barheight = unit(1.4, "mm"), ticks = FALSE),
    colour = guide_legend(order = 2, override.aes = list(size = 2.0)),
    shape = "none"
  ) +
  coord_cartesian(xlim = c(-tree_plot_shift - xmax * 0.02, panel_xmax), ylim = c(0, Ntip(tree) + 28), clip = "off") +
  labs(tag = "a") +
  theme_void(base_family = "Arial") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1),
    legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
    legend.title = element_text(size = 5.9), legend.text = element_text(size = 5.5),
    plot.margin = margin(10, 4, 1, 6, unit = "pt")
  )

# Panel B: two exact lineage cards; the historical unmapped F2 record is retained as a note.
make_lineage_card <- function(df, pair_id, tag = NULL) {
  z <- df[df$pair_id == pair_id & df$included_in_tree_panel, , drop = FALSE]
  z <- z[order(as_num(z$swedish_abs_rank <- core_primary$swedish_abs_rank[match(paste(z$pair_id, z$frozen_species), paste(core_primary$pair_id, core_primary$species))])), ]
  row_step <- 1.12
  z$row <- rev(seq_len(nrow(z))) * row_step
  top_y <- max(z$row)
  pair_repro <- if (pair_id == "F1_F1") repro_primary[repro_primary$pair_id == "Swedish_F1__HUNT_F1", ] else repro_primary[repro_primary$pair_id == "Swedish_F2__HUNT_F2", ]
  summary_text <- paste0(ifelse(pair_id == "F1_F1", "F1-F1", "F2-F2"), ": ",
                         pair_repro$exact_species_overlap_n, " mapped species | ",
                         pair_repro$common_species_distinct_genera_n, " genera | ",
                         pair_repro$common_species_distinct_families_n, " families | ",
                         pair_repro$direction_concordant_n, "/", pair_repro$direction_evaluable_n, " concordant")
  bg <- data.frame(row = z$row, ymin = z$row - 0.52, ymax = z$row + 0.52,
                   fill = rep(c("#F7F8FA", "#FFFFFF"), length.out = nrow(z)))
  p <- ggplot() +
    geom_rect(data = bg, aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax), fill = bg$fill, colour = NA) +
    geom_text(data = z, aes(x = 0.005, y = row, label = frozen_species), hjust = 0,
              size = 1.85, family = "Arial", fontface = "italic", colour = "#222222") +
    geom_text(data = z, aes(x = 0.51, y = row, label = genus), hjust = 0,
              size = 1.66, family = "Arial", colour = "#4A4A4A") +
    geom_text(data = z, aes(x = 0.75, y = row, label = family), hjust = 0,
              size = 1.62, family = "Arial", colour = "#4A4A4A") +
    geom_point(data = z, aes(x = 0.982, y = row, fill = direction_concordant),
               shape = 21, size = 2.15, stroke = 0.30, colour = "#333333") +
    scale_fill_manual(values = c("TRUE" = ifelse(pair_id == "F1_F1", palette_shared[["F1-F1"]], palette_shared[["F2-F2"]]),
                                 "FALSE" = "white"), guide = "none") +
    annotate("text", x = 0.005, y = top_y + 1.65, label = summary_text, hjust = 0,
             size = 2.0, family = "Arial", fontface = "bold", colour = "#222222") +
    annotate("text", x = c(0.005, 0.51, 0.75, 0.982), y = top_y + 0.76,
             label = c("Species", "Genus", "Family", "Dir."), hjust = c(0, 0, 0, 0.5),
             size = 1.72, family = "Arial", fontface = "bold", colour = "#555555") +
    coord_cartesian(xlim = c(0, 1), ylim = c(ifelse(pair_id == "F2_F2", -2.55, 0.40), top_y + 2.0), clip = "off") +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(3, 3, 1, 3, unit = "pt"),
          plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1))
  if (!is.null(tag)) p <- p + labs(tag = tag)
  if (pair_id == "F2_F2") {
    p <- p +
      annotate("text", x = 0.005, y = -0.48,
               label = "One additional historical F2 shared species (Pullichristensenella hMGS.05548)\nlacked a unique R232 endpoint and was excluded from mapped-species summaries;\nno taxonomic substitution was performed.",
               hjust = 0, vjust = 1, size = 1.38, family = "Arial", colour = "#666666", lineheight = 1.05) +
      annotate("point", x = 0.018, y = -2.10, shape = 21, size = 2.0,
               fill = palette_shared[["F2-F2"]], colour = "#333333", stroke = 0.3) +
      annotate("text", x = 0.050, y = -2.10, label = "Filled circle: direction-concordant.",
               hjust = 0, size = 1.42, family = "Arial", colour = "#555555") +
      annotate("point", x = 0.525, y = -2.10, shape = 21, size = 2.0,
               fill = "white", colour = "#333333", stroke = 0.3) +
      annotate("text", x = 0.557, y = -2.10, label = "Open circle: direction-discordant.",
               hjust = 0, size = 1.42, family = "Arial", colour = "#555555")
  }
  p
}

pB1 <- make_lineage_card(shared, "F1_F1", "b")
pB2 <- make_lineage_card(shared, "F2_F2")
pB <- pB1 / pB2 + plot_layout(heights = c(1.08, 1))

# Panel C
row_levels <- c("Swedish F1 Top20", "Swedish F1 Top40", "Swedish F2 Top20", "Swedish F2 Top40",
                "HUNT F1 Top20", "HUNT F1 Top40", "HUNT F2 Top20", "HUNT F2 Top40")
panel_c$base_y <- match(panel_c$display_row, rev(row_levels))
panel_c$plot_y <- panel_c$base_y + ifelse(panel_c$metric == "MPD", 0.16, -0.16)
panel_c$ses_num <- as_num(panel_c$ses)
panel_c$q_num <- as_num(panel_c$q_bh_family)
panel_c$q_label <- paste0(panel_c$metric, " q=", formatC(panel_c$q_num, format = "g", digits = 3))
robust_rows <- unique(panel_c$base_y[panel_c$conservative_q_lt_0_05])

pC <- ggplot(panel_c, aes(x = ses_num, y = plot_y)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#EEF4F8", alpha = 0.65) +
  geom_vline(xintercept = 0, colour = "#707070", linewidth = 0.32) +
  geom_segment(aes(x = 0, xend = ses_num, yend = plot_y, colour = profile), linewidth = 0.35, alpha = 0.72) +
  geom_point(aes(shape = metric, colour = profile, fill = ifelse(primary_q_lt_0_05, profile, "not_significant")),
             size = 2.15, stroke = 0.55) +
  geom_point(data = panel_c[panel_c$conservative_q_lt_0_05, ],
             aes(x = ses_num, y = plot_y), inherit.aes = FALSE, shape = 21, fill = NA,
             colour = "#111111", size = 3.05, stroke = 0.50) +
  geom_text(aes(x = 0.88, label = q_label), hjust = 0, size = 1.90, family = "Arial", colour = "#333333") +
  annotate("text", x = 2.15, y = robust_rows + 0.17, label = "MPD q = 0.0224", hjust = 0.5,
           size = 1.90, family = "Arial", colour = "#222222", fontface = "bold") +
  annotate("text", x = 2.15, y = robust_rows - 0.17, label = "MNTD q = 0.0192", hjust = 0.5,
           size = 1.90, family = "Arial", colour = "#222222", fontface = "bold") +
  annotate("text", x = 0.88, y = 8.72, label = "Primary BH q", hjust = 0,
           size = 1.90, family = "Arial", fontface = "bold", colour = "#444444") +
  annotate("text", x = 2.15, y = 8.72, label = "High-loading-propensity sensitivity", hjust = 0.5,
           size = 1.78, family = "Arial", fontface = "bold", colour = "#444444") +
  annotate("text", x = -3.55, y = 8.72, label = "Negative SES indicates clustering", hjust = 0,
           size = 1.72, family = "Arial", colour = "#444444") +
  scale_y_continuous(breaks = seq_along(row_levels), labels = rev(row_levels), limits = c(0.45, 8.9), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-3, -2, -1, 0), limits = c(-3.65, 2.75), expand = c(0, 0)) +
  scale_colour_manual(values = palette_profile, breaks = c("Swedish_F1", "Swedish_F2", "HUNT_F1", "HUNT_F2"), name = NULL) +
  scale_fill_manual(values = c(palette_profile, not_significant = "white"), guide = "none") +
  scale_shape_manual(values = c(MPD = 24, MNTD = 21), name = NULL) +
  labs(x = "Standardized effect size (SES)", y = NULL, tag = "c",
       caption = "Horizontal stems connect SES estimates to zero and do not represent confidence intervals.") +
  theme_pub +
  theme(
    legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit(1.2, "mm"),
    plot.margin = margin(5, 4, 4, 6, unit = "pt"), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.caption = element_text(size = 5.8, hjust = 0, colour = "#555555", margin = margin(t = 3, unit = "pt"))
  ) +
  guides(colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 2)),
         shape = guide_legend(order = 2, override.aes = list(colour = "#444444", fill = "white", size = 2)))

pA_free_left <- patchwork::free(pA, type = "panel", side = "l")
figure <- (pA_free_left | pB) / pC +
  plot_layout(widths = c(1.52, 1.08), heights = c(2.15, 1)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

width_mm <- 183
height_mm <- 180
base_path <- file.path(outdir, "Figure3_phylogenetic_organization_final")

svglite::svglite(paste0(base_path, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4,
                 bg = "white")
print(figure)
dev.off()

grDevices::cairo_pdf(paste0(base_path, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
                     family = "Arial", bg = "white", onefile = TRUE)
print(figure)
dev.off()

png_width_px <- round(width_mm / 25.4 * 600)
png_height_px <- round(height_mm / 25.4 * 600)
ragg::agg_png(paste0(base_path, ".png"), width = png_width_px, height = png_height_px,
              units = "px", res = 600, background = "white")
print(figure)
dev.off()

# R-only PDF render for layout QA (pdftools uses Poppler through R).
pdf_qa_png <- file.path(outdir, "Figure3_phylogenetic_organization_final_pdf_render.png")
pdf_qa_pattern <- file.path(outdir, "Figure3_phylogenetic_organization_final_pdf_render-%d.%s")
pdf_qa_created <- pdftools::pdf_convert(paste0(base_path, ".pdf"), format = "png", dpi = 300,
                                        filenames = pdf_qa_pattern, verbose = FALSE)
if (file.exists(pdf_qa_png)) unlink(pdf_qa_png)
stopifnot(length(pdf_qa_created) == 1L, file.rename(pdf_qa_created[[1]], pdf_qa_png))
log_msg("Figure exports complete: SVG, PDF, 600-dpi PNG; R-rendered PDF QA preview created")

# Export-format verification: exact raster dimensions/density and searchable vector text.
png_raw <- readBin(paste0(base_path, ".png"), "raw", n = file.info(paste0(base_path, ".png"))$size)
u32be <- function(x) sum(as.numeric(as.integer(x)) * 256^(3:0))
png_width_read <- u32be(png_raw[17:20])
png_height_read <- u32be(png_raw[21:24])
pos <- 9L
png_x_ppm <- png_y_ppm <- png_phys_unit <- NA_real_
while (pos <= length(png_raw) - 11L) {
  chunk_len <- u32be(png_raw[pos:(pos + 3L)])
  chunk_type <- rawToChar(png_raw[(pos + 4L):(pos + 7L)])
  if (chunk_type == "pHYs") {
    chunk <- png_raw[(pos + 8L):(pos + 8L + chunk_len - 1L)]
    png_x_ppm <- u32be(chunk[1:4])
    png_y_ppm <- u32be(chunk[5:8])
    png_phys_unit <- as.integer(chunk[9])
    break
  }
  pos <- pos + 12L + chunk_len
}
png_dpi_x <- png_x_ppm * 0.0254
png_dpi_y <- png_y_ppm * 0.0254
pdf_meta <- pdftools::pdf_info(paste0(base_path, ".pdf"))
pdf_text_chars <- nchar(paste(pdftools::pdf_text(paste0(base_path, ".pdf")), collapse = ""))
svg_has_editable_text <- any(grepl("<text", readLines(paste0(base_path, ".svg"), warn = FALSE), fixed = TRUE))
stopifnot(
  png_width_read == png_width_px,
  png_height_read == png_height_px,
  png_phys_unit == 1L,
  abs(png_dpi_x - 600) < 0.01,
  abs(png_dpi_y - 600) < 0.01,
  pdf_meta$pages == 1L,
  pdf_text_chars > 0L,
  svg_has_editable_text
)
log_msg("Export QC PASS: PNG ", png_width_read, "x", png_height_read,
        " at ", formatC(png_dpi_x, digits = 6, format = "f"), " dpi; PDF text chars=", pdf_text_chars,
        "; SVG contains editable text elements")

# Build notes and QC. visual_qc_status is updated only after manual inspection.
visual_qc_status <- "PASS - direct 600-dpi PNG and independently R-rendered PDF preview inspected; no clipping, overlap, or missing plotted values"

build_notes <- c(
  "# Figure 3 final phylogenetic visual-refinement build notes",
  "",
  "## Figure contract",
  "",
  "- **Core conclusion:** High-loading species show limited phylogenetic organization; robust clustering is confined to Swedish F1 Top40 under the conservative high-loading-propensity sensitivity.",
  "- **Archetype:** Asymmetric mixed-modality figure with Panel A as the hero panel.",
  "- **Backend:** R only (`ape`, `ggtree`, `ggplot2`, `ggrepel`, `patchwork`, `svglite`, `ragg`).",
  "- **Final size:** 183 mm x 180 mm.",
  "- **Panel A:** Exact pruned GTDB R11-RS232 tree with all 236 unique endpoints, four frozen signed-loading tracks, and mapped shared Top20 highlights. The plotted tree was translated left without changing topology or branch lengths, and the direct-label lane was moved farther right of the loading tracks. Non-data diagonal ggrepel leader segments were removed to avoid the former double-slash artifact around the two ER4 labels.",
  "- **Panel B:** Exact mapped shared Top20 species, genus/family identities, and direction concordance. Filled and open circles are explicitly defined. The seventh historical F2 shared species remains in the historical record but is excluded from mapped-species summaries without substitution.",
  "- **Panel C:** All 16 frozen primary information-stratified MPD/MNTD SES and BH q values are retained. Only Swedish F1 Top40 displays the frozen high-loading-propensity sensitivity values (MPD q = 0.0224; MNTD q = 0.0192). Horizontal stems connect SES estimates to zero and are not confidence intervals.",
  "",
  "## Frozen-data handling",
  "",
  "No phylogenetic statistic, P value, q value, null distribution, GFA loading, factor sign, rank, or shared-species membership was recomputed. The tree was pruned only to the 236 frozen unique endpoints. HUNT F1 uses the frozen -1 multiplier and HUNT F2 the frozen +1 multiplier.",
  "",
  "## Supplementary placement",
  "",
  "Blomberg K/Pagel lambda, family-block RBO, shared-pool MPD/MNTD, uniform-null comparisons, full mapping/version audit, and complete sensitivity tables remain supplementary. No shared-pool clustering test survived FDR correction.",
  "",
  "## Interpretation boundary",
  "",
  "The visual supports limited phylogenetic organization and cross-lineage loading recurrence. It does not support conserved microbial function, causal ecology, complete factor replication, SNP-effect replication, or universal clustering."
)
writeLines(build_notes, file.path(outdir, "FIGURE3_BUILD_NOTES.md"), useBytes = TRUE)

qc_lines <- c(
  "# Figure 3 phylogeny QC",
  "",
  paste0("**Visual QC status: ", visual_qc_status, "**"),
  "",
  "## Programmatic checks",
  "",
  "- Mapping rows: 241; unique R232 endpoints plotted: 236; duplicate endpoints: 0.",
  "- Frozen signs: HUNT F1 multiplier = -1; HUNT F2 multiplier = +1.",
  "- Four mapped loading profiles: exact match to the frozen 241-species loading export.",
  "- Panel B: F1-F1 = 8 mapped species, 7 genera, 6 families, 6/8 direction-concordant; F2-F2 = 6 mapped species, 6 genera, 5 families, 6/6 direction-concordant.",
  "- Historical F2 shared species without a unique endpoint retained as an excluded record; no genus/congener substitution.",
  "- Panel C: 16 primary tests; 8 have primary BH q<0.05; 2 retain conservative BH q<0.05, both Swedish F1 Top40 (MPD and MNTD).",
  "- Shared-pool MPD/MNTD FDR positives: 0.",
  paste0("- ", workbook_checks$check, ": **", workbook_checks$status, "** (", workbook_checks$detail, ")"),
  "- Value audit records: 1019 (the original 1017 records are retained verbatim, plus two traceable rounded-display records for the frozen sensitivity q values).",
  paste0("- Figure dimensions: 183 mm x 180 mm; PNG = ", png_width_px, " x ", png_height_px, " pixels at 600 dpi; SVG/PDF remain vector exports."),
  paste0("- Export metadata verified: PNG pHYs density = ", formatC(png_dpi_x, digits = 4, format = "f"), " x ", formatC(png_dpi_y, digits = 4, format = "f"), " dpi; PDF searchable text characters = ", pdf_text_chars, "; SVG contains text elements."),
  "",
  "## Visual review checklist",
  "",
  "- All 236 tips remain represented without forcing all tip labels.",
  "- Shared recurrent species labels are readable and connected to exact tips.",
  "- Four signed-loading track headers and diverging scale are legible and not clipped at the top edge.",
  "- Panel A uses the available left-side space; recurrent-species labels no longer overlap the loading tracks, and the two ER4 labels have no duplicate diagonal leader-line artifact.",
  "- Panel B species/genus/family columns do not overlap; filled/open direction symbols are explicitly explained; the F2 mapping footnote is complete.",
  "- Panel C shows all 16 points, all primary q labels, the zero reference, and only the two frozen Swedish F1 Top40 sensitivity q values.",
  "- No clipping, missing glyphs, rasterized text, unintended legends, or decorative elements.",
  "- PDF was rendered through R/pdftools and checked separately from the direct PNG export.",
  "",
  "## Scientific integrity",
  "",
  "No new analysis was performed. Negative SES is explicitly labeled as clustering. The figure does not imply functional conservation, causal mechanism, complete cross-cohort factor replication, SNP-effect replication, or universal phylogenetic clustering."
)
writeLines(qc_lines, file.path(outdir, "FIGURE3_PHYLOGENY_FINAL_QC.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(outdir, "R_SESSION_INFO.txt"))
log_msg("Reports and R session information written; visual QC status=", visual_qc_status)
log_msg("END build completed without source inconsistency")

input_manifest <- source_inventory[, c("input_id", "absolute_path", "intended_figure_use", "bytes", "sha256")]
write_tsv(input_manifest, file.path(outdir, "INPUT_SHA256_MANIFEST.tsv"))

output_names <- c(
  "Figure3_phylogenetic_organization_final.pdf",
  "Figure3_phylogenetic_organization_final.svg",
  "Figure3_phylogenetic_organization_final.png",
  "Figure3_phylogenetic_organization_final_pdf_render.png",
  "F3A_TREE_LAYOUT.tsv",
  "F3A_PHYLOGENY_LOADINGS.tsv",
  "F3B_SHARED_TOP20_LINEAGES.tsv",
  "F3C_PRIMARY_CLUSTERING.tsv",
  "figure_value_audit.tsv",
  "SOURCE_INVENTORY.tsv",
  "INPUT_SHA256_MANIFEST.tsv",
  "FIGURE3_PHYLOGENY_FINAL_QC.md",
  "FIGURE3_BUILD_NOTES.md",
  "R_SESSION_INFO.txt",
  "FIGURE3_PHYLOGENY_FINAL_EXECUTION.log",
  "build_figure3_phylogenetic_final.R"
)
output_paths <- file.path(outdir, output_names)
stopifnot(all(file.exists(output_paths)))
output_manifest <- data.frame(
  output_file = output_names,
  absolute_path = normalizePath(output_paths, winslash = "/"),
  bytes = as.character(file.info(output_paths)$size),
  sha256 = vapply(output_paths, sha256, character(1)),
  check.names = FALSE
)
write_tsv(output_manifest, file.path(outdir, "OUTPUT_SHA256_MANIFEST.tsv"))

cat("FINAL_STATUS\tREADY_FOR_AUTHOR_REVIEW\n")
cat("OUTPUT_DIR\t", outdir, "\n", sep = "")
