if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
root <- file.path(Sys.getenv("MHFNDD_PROJECT_ROOT"), "FULL_SCALE_GFA_FINAL_EVALUATION")
model_dir <- file.path(root, "models")
out_dir <- file.path(root, "cross_cohort_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(file.path(root, "inputs_241_manifest.tsv"), check.names = FALSE,
                       stringsAsFactors = FALSE)
snp_manifest <- read.delim(file.path(root, "common_genomewide_ld_controlled_training_manifest.tsv"),
                           check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 241L, nrow(snp_manifest) == 330704L)

extract_one <- function(cohort) {
  path <- file.path(model_dir, sprintf("%s241_genomewide_ldcontrolled_Rldsc_seed20261014.rds", cohort))
  x <- readRDS(path)
  stopifnot(identical(dim(x$F_hat), c(241L, 2L)),
            identical(dim(x$L_hat), c(330704L, 2L)),
            all(is.finite(x$F_hat)), all(is.finite(x$L_hat)))
  rownames(x$F_hat) <- manifest$swedish_species
  rownames(x$L_hat) <- snp_manifest$canonical_key
  small <- list(
    cohort = cohort,
    method = x$method,
    F_hat = x$F_hat,
    L_hat = x$L_hat,
    F_hat_single = x$F_hat_single,
    F_hat_multi = x$F_hat_multi,
    scale = x$scale,
    gfa_pve = x$gfa_pve,
    params = x$params,
    R = x$R
  )
  saveRDS(small, file.path(out_dir, sprintf("%s241_Rldsc_primary_factor_core.rds", cohort)),
          compress = "xz")
  sizes <- vapply(small, object.size, numeric(1))
  write.table(data.frame(field = names(sizes), bytes = as.numeric(sizes)),
              file.path(out_dir, sprintf("%s241_core_field_sizes.tsv", cohort)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines(capture.output(str(small, max.level = 2L, give.attr = FALSE)),
             file.path(out_dir, sprintf("%s241_core_structure.txt", cohort)))
  rm(x, small); invisible(gc())
}

extract_one("Swedish")
extract_one("HUNT")

sw <- readRDS(file.path(out_dir, "Swedish241_Rldsc_primary_factor_core.rds"))
hu <- readRDS(file.path(out_dir, "HUNT241_Rldsc_primary_factor_core.rds"))

pearson_F <- cor(sw$F_hat, hu$F_hat, method = "pearson")
spearman_F <- cor(sw$F_hat, hu$F_hat, method = "spearman")
pearson_L <- cor(sw$L_hat, hu$L_hat, method = "pearson")
spearman_L <- cor(sw$L_hat, hu$L_hat, method = "spearman")

# With two factors, enumerate both possible one-to-one assignments and choose the
# assignment maximizing the summed absolute species-loading correlation.
perms <- rbind(c(1L, 2L), c(2L, 1L))
scores <- apply(perms, 1L, function(p) sum(abs(pearson_F[cbind(1:2, p)])))
best <- perms[which.max(scores), ]
signs <- sign(pearson_F[cbind(1:2, best)])
signs[signs == 0] <- 1

pairing <- data.frame(
  swedish_factor = 1:2,
  hunt_factor = best,
  sign_applied_to_hunt = signs,
  species_pearson = pearson_F[cbind(1:2, best)],
  species_spearman = spearman_F[cbind(1:2, best)],
  snp_pearson = pearson_L[cbind(1:2, best)],
  snp_spearman = spearman_L[cbind(1:2, best)],
  stringsAsFactors = FALSE
)
pairing$species_pearson_aligned <- pairing$species_pearson * signs
pairing$species_spearman_aligned <- pairing$species_spearman * signs
pairing$snp_pearson_aligned <- pairing$snp_pearson * signs
pairing$snp_spearman_aligned <- pairing$snp_spearman * signs
write.table(pairing, file.path(out_dir, "primary_factor_optimal_pairing.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

matrix_to_long <- function(m, metric, layer) {
  z <- expand.grid(swedish_factor = seq_len(nrow(m)),
                   hunt_factor = seq_len(ncol(m)))
  z$value <- as.vector(m)
  z$metric <- metric; z$layer <- layer
  z[, c("layer", "metric", "swedish_factor", "hunt_factor", "value")]
}
all_corr <- rbind(matrix_to_long(pearson_F, "pearson", "species_loadings"),
                  matrix_to_long(spearman_F, "spearman", "species_loadings"),
                  matrix_to_long(pearson_L, "pearson", "training_snp_loadings"),
                  matrix_to_long(spearman_L, "spearman", "training_snp_loadings"))
write.table(all_corr, file.path(out_dir, "all_factor_cross_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

factor_summary <- function(core, cohort) {
  do.call(rbind, lapply(seq_len(ncol(core$F_hat)), function(k) {
    f <- core$F_hat[, k]; l <- core$L_hat[, k]
    data.frame(cohort = cohort, factor = k,
      species_nonzero_exact = sum(f != 0), species_abs_gt_1e8 = sum(abs(f) > 1e-8),
      species_l2_norm = sqrt(sum(f^2)), species_abs_median = median(abs(f)),
      species_abs_max = max(abs(f)), snp_nonzero_exact = sum(l != 0),
      snp_abs_gt_1e8 = sum(abs(l) > 1e-8), snp_l2_norm = sqrt(sum(l^2)),
      snp_abs_median = median(abs(l)), snp_abs_max = max(abs(l)))
  }))
}
summ <- rbind(factor_summary(sw, "Swedish"), factor_summary(hu, "HUNT"))
write.table(summ, file.path(out_dir, "primary_factor_numeric_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

aligned_species <- data.frame(trait_order = manifest$trait_order,
                              species = manifest$swedish_species,
                              stringsAsFactors = FALSE)
aligned_snps <- snp_manifest[, c("chrom", "position_hg38", "allele1", "allele2",
                                 "reference_snp_id", "canonical_key")]
for (k in 1:2) {
  j <- best[k]; s <- signs[k]
  aligned_species[[sprintf("Swedish_F%d", k)]] <- sw$F_hat[, k]
  aligned_species[[sprintf("HUNT_F%d_aligned", j)]] <- hu$F_hat[, j] * s
  aligned_snps[[sprintf("Swedish_L%d", k)]] <- sw$L_hat[, k]
  aligned_snps[[sprintf("HUNT_L%d_aligned", j)]] <- hu$L_hat[, j] * s
}
write.table(aligned_species, file.path(out_dir, "aligned_species_loadings.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(aligned_snps, file.path(out_dir, "aligned_training_snp_loadings.tsv.gz"),
            sep = "\t", row.names = FALSE, quote = FALSE)

top_rows <- do.call(rbind, lapply(1:2, function(k) {
  j <- best[k]; s <- signs[k]
  ix <- order(abs(sw$F_hat[, k]), decreasing = TRUE)[1:25]
  data.frame(swedish_factor = k, hunt_factor = j, rank_by_swedish_abs = 1:25,
             trait_order = manifest$trait_order[ix], species = manifest$swedish_species[ix],
             swedish_loading = sw$F_hat[ix, k], hunt_loading_aligned = hu$F_hat[ix, j] * s,
             stringsAsFactors = FALSE)
}))
write.table(top_rows, file.path(out_dir, "top25_species_by_swedish_factor.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(c(
  "IMPORTANT_INTERPRETATION_LIMITATION",
  "Species-loading correlations assess structural replication only.",
  "Training-SNP loading correlations reuse the same frozen SNP set and are descriptive, not held-out validation.",
  "Neither quantity alone satisfies Level 2 CROSS_COHORT_FACTOR_SUPPORTED.",
  "A pre-defined LD-independent held-out genetic evaluation or calibrated SNP statistic remains required."
), file.path(out_dir, "INTERPRETATION_LIMITATIONS.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "primary_factor_comparison_sessionInfo.txt"))
