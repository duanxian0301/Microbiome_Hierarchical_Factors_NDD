if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
suppressPackageStartupMessages(library(coloc))

root <- file.path(Sys.getenv("MHFNDD_PROJECT_ROOT"), "PHASE16_SCZ_THREE_WAY_SIGNAL_AUDIT")
harm_path <- file.path(root, "THREE_WAY_VARIANT_HARMONIZATION.tsv")
common_path <- file.path(root, "COMMON_SNP_UNIVERSE.tsv")
phase14_path <- file.path(Sys.getenv("MHFNDD_PROJECT_ROOT"), "PHASE14_TARGETED_VALIDATION", "B_SCZ_CD_UC_LOCI", "coloc_inputs", "SCZ_3_17582158_18082158__Pullichristensenella_excrementipullorum__coloc_input.tsv")
harm <- read.delim(harm_path, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("NA", ""))
common <- read.delim(common_path, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("NA", ""))
phase14 <- read.delim(phase14_path, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("NA", ""))

p12_values <- c(1e-6, 1e-5, 1e-4)
summary_rows <- list()
snp_rows <- list()
warnings_seen <- character()

make_scz <- function(x, include_n=FALSE) {
  ncase <- median(x$scz_ncase, na.rm=TRUE)
  ncontrol <- median(x$scz_ncontrol, na.rm=TRUE)
  d <- list(
    beta=x$scz_harmonized_beta,
    varbeta=x$scz_se^2,
    snp=x$snp,
    position=x$position_grch38,
    type="cc",
    s=ncase/(ncase+ncontrol)
  )
  if (include_n) d$N <- ncase+ncontrol
  d
}

make_species <- function(x) {
  list(
    beta=x$species_harmonized_beta,
    varbeta=x$species_se^2,
    snp=x$snp,
    position=x$position_grch38,
    type="quant",
    N=16017,
    MAF=x$species_reference_maf_1000g_eur
  )
}

make_eqtl <- function(x) {
  list(
    beta=x$eqtl_harmonized_beta,
    varbeta=x$eqtl_se_reconstructed^2,
    snp=x$snp,
    position=x$position_grch38,
    type="quant",
    sdY=1
  )
}

run_fit <- function(x, cell_type, pair, scope, d1, d2, input_quality) {
  stopifnot(length(d1$beta) == nrow(x), length(d2$beta) == nrow(x))
  stopifnot(!anyDuplicated(x$snp))
  for (p12 in p12_values) {
    fit <- tryCatch(
      withCallingHandlers(
        coloc.abf(d1, d2, p1=1e-4, p2=1e-4, p12=p12),
        warning=function(w) {
          warnings_seen <<- c(warnings_seen, paste(cell_type, pair, scope, format(p12, scientific=TRUE), conditionMessage(w), sep="\t"))
          invokeRestart("muffleWarning")
        }
      ),
      error=function(e) e
    )
    if (inherits(fit, "error")) {
      summary_rows[[length(summary_rows)+1]] <<- data.frame(
        cell_type=cell_type, comparison=pair, analysis_scope=scope,
        p1=1e-4, p2=1e-4, p12=p12, snp_n=nrow(x), status="ERROR",
        PP_H0=NA, PP_H1=NA, PP_H2=NA, PP_H3=NA, PP_H4=NA,
        top_snp=NA, top_position_grch38=NA, top_snp_PP_H4=NA,
        rs2337387_SNP_PP_H4=NA, rs7428585_species_lead_SNP_PP_H4=NA,
        input_quality=input_quality, error=conditionMessage(fit), stringsAsFactors=FALSE
      )
      next
    }
    s <- fit$summary
    result <- fit$results
    top_i <- which.max(result$SNP.PP.H4)
    rs_i <- which(result$snp == "rs2337387")
    species_i <- which(result$snp == "rs7428585")
    summary_rows[[length(summary_rows)+1]] <<- data.frame(
      cell_type=cell_type, comparison=pair, analysis_scope=scope,
      p1=1e-4, p2=1e-4, p12=p12, snp_n=as.integer(s[["nsnps"]]), status="COMPLETE",
      PP_H0=as.numeric(s[["PP.H0.abf"]]), PP_H1=as.numeric(s[["PP.H1.abf"]]),
      PP_H2=as.numeric(s[["PP.H2.abf"]]), PP_H3=as.numeric(s[["PP.H3.abf"]]),
      PP_H4=as.numeric(s[["PP.H4.abf"]]),
      top_snp=as.character(result$snp[[top_i]]), top_position_grch38=as.integer(result$position[[top_i]]),
      top_snp_PP_H4=as.numeric(result$SNP.PP.H4[[top_i]]),
      rs2337387_SNP_PP_H4=if(length(rs_i)==1) as.numeric(result$SNP.PP.H4[[rs_i]]) else NA,
      rs7428585_species_lead_SNP_PP_H4=if(length(species_i)==1) as.numeric(result$SNP.PP.H4[[species_i]]) else NA,
      input_quality=input_quality, error="", stringsAsFactors=FALSE
    )
    snp_rows[[length(snp_rows)+1]] <<- data.frame(
      cell_type=cell_type, comparison=pair, analysis_scope=scope,
      p1=1e-4, p2=1e-4, p12=p12,
      snp=as.character(result$snp), position_grch38=as.integer(result$position),
      SNP_PP_H4=as.numeric(result$SNP.PP.H4),
      is_rs2337387=as.character(result$snp)=="rs2337387",
      is_species_lead=as.character(result$snp)=="rs7428585",
      rank_SNP_PP_H4=rank(-as.numeric(result$SNP.PP.H4), ties.method="min"),
      stringsAsFactors=FALSE
    )
  }
}

# Exact reproduction of the historical Phase14 SCZ-species setup on 1,030 SNPs.
x14 <- data.frame(
  snp=phase14$SNP_KEY,
  position_grch38=phase14$BP,
  species_harmonized_beta=phase14$SPECIES_BETA_ALIGNED,
  species_se=phase14$SPECIES_SE,
  species_reference_maf_1000g_eur=phase14$REFERENCE_MAF_1000G_EUR,
  scz_harmonized_beta=phase14$DISEASE_BETA,
  scz_se=phase14$DISEASE_SE,
  scz_ncase=phase14$DISEASE_NCASE,
  scz_ncontrol=phase14$DISEASE_NCONTROL,
  stringsAsFactors=FALSE
)
run_fit(
  x14, "not_applicable", "SCZ_vs_species", "phase14_reproduction_1030",
  make_species(x14), make_scz(x14, include_n=TRUE),
  "Historical coverage reproduction; reference MAF and study-level species N"
)

cells <- c("Excitatory.neurons", "OPCs...COPs")
for (cell in cells) {
  x_common <- common[common$cell_type == cell, ]
  x_612 <- harm[harm$cell_type == cell & (is.na(harm$exclusion_reason) | harm$exclusion_reason == "palindromic_without_source_specific_AF"), ]
  x_612 <- x_612[is.finite(x_612$scz_harmonized_beta) & is.finite(x_612$species_harmonized_beta) & is.finite(x_612$eqtl_harmonized_beta), ]

  # Primary, strictly harmonized three-way universe.
  run_fit(
    x_common, cell, "SCZ_vs_species", "three_way_common_501",
    make_species(x_common), make_scz(x_common, include_n=TRUE),
    "Primary common universe; non-palindromic; species reference MAF/study-level N sensitivity"
  )
  run_fit(
    x_common, cell, "SCZ_vs_TBC1D5_eQTL", "three_way_common_501",
    make_scz(x_common, include_n=FALSE), make_eqtl(x_common),
    "Primary common universe; eQTL SE reconstructed from beta/P; eQTL sdY=1"
  )
  run_fit(
    x_common, cell, "species_vs_TBC1D5_eQTL", "three_way_common_501",
    make_species(x_common), make_eqtl(x_common),
    "Primary common universe; both species and eQTL metadata limitations retained"
  )

  # Historical-coverage QC: mirrors Phase15's 612 overlapping SNPs, including
  # palindromic variants whose source-specific allele frequencies are absent.
  run_fit(
    x_612, cell, "SCZ_vs_TBC1D5_eQTL", "phase15_coverage_reproduction_612",
    make_scz(x_612, include_n=FALSE), make_eqtl(x_612),
    "Coverage QC only; includes palindromic variants and is not the Phase16 primary result"
  )
  run_fit(
    x_612, cell, "species_vs_TBC1D5_eQTL", "phase15_coverage_reproduction_612",
    make_species(x_612), make_eqtl(x_612),
    "Coverage QC only; includes palindromic variants and is not the Phase16 primary result"
  )
}

summary_out <- do.call(rbind, summary_rows)
snp_out <- do.call(rbind, snp_rows)
write.table(summary_out, file.path(root, "PAIRWISE_COLOC_REANALYSIS.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="NA")
write.table(snp_out, file.path(root, "SNP_LEVEL_POSTERIOR_COMPARISON.tsv"), sep="\t", row.names=FALSE, quote=FALSE, na="NA")
writeLines(unique(warnings_seen), file.path(root, "logs", "phase16_coloc_warnings.log"))
cat(sprintf("summary_rows=%d snp_rows=%d complete=%d\n", nrow(summary_out), nrow(snp_out), sum(summary_out$status=="COMPLETE")))
