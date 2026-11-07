if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(data.table); library(digest)})

ROOT <- Sys.getenv("MHFNDD_PROJECT_ROOT")
BASE <- file.path(ROOT, "FULL_SCALE_GFA_FINAL_EVALUATION")
P9 <- file.path(BASE, "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA")
P7 <- file.path(BASE, "PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY")
B <- 20000L
SEED <- 20261024L
RBO_P <- 0.95

rd <- function(path) fread(path, sep = "\t", data.table = FALSE, check.names = FALSE, na.strings = c("NA", ""))
wr <- function(x, name) fwrite(x, file.path(P9, name), sep = "\t", quote = FALSE, na = "NA", logical01 = FALSE)
cosine <- function(a, b) {
  den <- sqrt(sum(a^2) * sum(b^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  sum(a * b) / den
}
rank_order <- function(v, trait_order, species) order(-abs(v), trait_order, species, method = "radix")
rbo_orders <- function(a, b, p = RBO_P) {
  k <- length(a); pos_a <- integer(k); pos_b <- integer(k)
  pos_a[a] <- seq_len(k); pos_b[b] <- seq_len(k)
  m <- pmax(pos_a, pos_b)
  w <- (1 - p) * p^(seq_len(k) - 1L) / seq_len(k)
  tail_sum <- rev(cumsum(rev(w)))
  sum(tail_sum[m]) + p^k
}
top_overlap <- function(a, b, k) length(intersect(head(a, k), head(b, k))) / k
direction_concordance <- function(a, b, order_a, order_b, k) {
  common <- intersect(head(order_a, k), head(order_b, k))
  if (!length(common)) return(NA_real_)
  mean(sign(a[common]) == sign(b[common]))
}

global <- rd(file.path(BASE, "BIDIRECTIONAL_SPECIES_LOADING_H2.tsv"))
species <- as.character(global$species)
trait_order <- as.integer(global$trait_order)
stopifnot(length(species) == 241L, !anyDuplicated(species))
G_SW <- cbind(F1 = global$Swedish_F1_loading_raw, F2 = global$Swedish_F2_loading_raw)
G_HU <- cbind(F1 = global$HUNT_loading_aligned_to_Swedish_F1, F2 = global$HUNT_loading_aligned_to_Swedish_F2)

stageb <- rd(file.path(P9, "STAGEB_HUNT_VALIDATION_MANIFEST.tsv"))
explore <- rd(file.path(P9, "STAGEB_POST_STAGEA_EXPLORATORY_MANIFEST.tsv"))
primary_ids <- as.character(stageb$annotation_id)
explore_ids <- as.character(explore$annotation_id)
stopifnot(length(primary_ids) == 18L, length(explore_ids) == 4L, !anyDuplicated(c(primary_ids, explore_ids)))
all_ids <- c(primary_ids, explore_ids)

ledger <- rd(file.path(P9, "PHASE9_BATCH_EXECUTION_LEDGER.tsv"))
sel_ledger <- ledger[ledger$analysis_set == "target" & ledger$annotation_id %in% all_ids & ledger$cohort %in% c("Swedish", "HUNT"), ]
stopifnot(nrow(sel_ledger) == 44L, all(sel_ledger$status == "SUCCESS"))

loadings <- list(Swedish = list(), HUNT = list())
pves <- list(Swedish = list(), HUNT = list())
summary_rows <- list()
sr <- 0L
for (cohort in c("Swedish", "HUNT")) {
  for (aid in all_ids) {
    task <- sel_ledger[sel_ledger$cohort == cohort & sel_ledger$annotation_id == aid, , drop = FALSE]
    stopifnot(nrow(task) == 1L, file.exists(task$summary_path))
    s <- rd(task$summary_path)
    x <- rd(s$loadings_path[1]); pv <- rd(s$pve_path[1])
    nf <- as.integer(s$retained_factors[1]); fc <- if (nf) paste0("F", seq_len(nf)) else character()
    stopifnot(nrow(x) == 241L, identical(as.character(x$species), species), nrow(pv) == 241L,
              identical(as.character(pv$species), species))
    L <- if (nf) as.matrix(x[, fc, drop = FALSE]) else matrix(numeric(), nrow = 241L, ncol = 0L)
    P <- if (nf) as.matrix(pv[, fc, drop = FALSE]) else matrix(numeric(), nrow = 241L, ncol = 0L)
    stopifnot(nf == ncol(L), nf == ncol(P), nf == 0L || all(is.finite(L)), nf == 0L || all(is.finite(P)))
    loadings[[cohort]][[aid]] <- L; pves[[cohort]][[aid]] <- P
    sr <- sr + 1L
    summary_rows[[sr]] <- data.frame(
      annotation_id = aid, analysis_family = if (aid %in% primary_ids) "PRIMARY_18" else "POST_STAGEA_EXPLORATORY",
      organ = ifelse(startsWith(aid, "gut__"), "gut", "brain"), cohort = cohort,
      model_status = s$model_status[1], variants = as.integer(s$variants[1]), traits = as.integer(s$traits[1]),
      sample_size = as.integer(s$sample_size[1]), retained_factors = nf,
      retained_multitrait_factors = as.integer(s$retained_multitrait_factors[1]),
      retained_singletrait_factors = as.integer(s$retained_singletrait_factors[1]),
      maxiter_reached = as.logical(s$maxiter_reached[1]), convergence_ok = !as.logical(s$maxiter_reached[1]),
      elapsed_seconds = as.numeric(s$elapsed_seconds[1]), input_z_path = s$input_z_path[1],
      input_z_sha256 = s$input_z_sha256[1],
      input_hash_matches_ledger_summary = tolower(s$input_z_sha256[1]) == tolower(task$input_z_sha256[1]),
      R_version = s$R_version[1], GFA_version = s$GFA_version[1], flashier_version = s$flashier_version[1], ebnm_version = s$ebnm_version[1],
      model_path = s$model_path[1], model_bytes = as.numeric(s$model_bytes[1]), model_sha256 = s$model_sha256[1],
      loadings_path = s$loadings_path[1], loadings_sha256 = digest(s$loadings_path[1], algo = "sha256", file = TRUE, serialize = FALSE),
      pve_path = s$pve_path[1], pve_sha256 = digest(s$pve_path[1], algo = "sha256", file = TRUE, serialize = FALSE),
      species_order_241_ok = TRUE, loadings_finite = nf == 0L || all(is.finite(L)), pve_finite = nf == 0L || all(is.finite(P)),
      completed_output_files_present = task$status[1] == "SUCCESS" && all(file.exists(c(s$model_path[1], s$loadings_path[1], s$pve_path[1]))),
      stringsAsFactors = FALSE
    )
  }
}
model_summary <- rbindlist(summary_rows, fill = TRUE)

match_cell <- function(L, G, cohort, aid, P) {
  nf <- ncol(L)
  factor_audit <- list(); metrics <- vector("list", nf)
  if (nf) for (j in seq_len(nf)) {
    l <- L[, j]
    c1 <- cosine(l, G[, 1]); c2 <- cosine(l, G[, 2])
    metrics[[j]] <- list(c1 = c1, c2 = c2, norm = sqrt(sum(l^2)),
                         pve_sum = sum(P[, j]), pve_median = median(P[, j]), pve_max = max(P[, j]))
    for (d in 1:2) {
      gd <- G[, d]; ord_l <- rank_order(l, trait_order, species); ord_g <- rank_order(gd, trait_order, species)
      factor_audit[[length(factor_audit) + 1L]] <- data.frame(
        annotation_id = aid, cohort = cohort, factor = paste0("F", j), global_direction = paste0("F", d),
        signed_cosine = if (d == 1) c1 else c2, absolute_cosine = abs(if (d == 1) c1 else c2),
        pearson = cor(l, gd), spearman = cor(l, gd, method = "spearman"),
        rbo_p095 = rbo_orders(ord_l, ord_g), top20_overlap = top_overlap(ord_l, ord_g, 20L),
        top40_overlap = top_overlap(ord_l, ord_g, 40L), loading_l2_norm = metrics[[j]]$norm,
        factor_pve_sum = metrics[[j]]$pve_sum, factor_pve_median = metrics[[j]]$pve_median,
        factor_pve_max = metrics[[j]]$pve_max, stringsAsFactors = FALSE)
    }
  }
  j1 <- j2 <- NA_integer_
  if (nf >= 2L) {
    cand <- expand.grid(j1 = seq_len(nf), j2 = seq_len(nf)); cand <- cand[cand$j1 != cand$j2, ]
    cand$sum <- mapply(function(x, y) abs(metrics[[x]]$c1) + abs(metrics[[y]]$c2), cand$j1, cand$j2)
    cand$minimum <- mapply(function(x, y) min(abs(metrics[[x]]$c1), abs(metrics[[y]]$c2)), cand$j1, cand$j2)
    cand <- cand[order(-cand$sum, -cand$minimum, cand$j1, cand$j2), ]
    j1 <- cand$j1[1]; j2 <- cand$j2[1]
  } else if (nf == 1L) {
    if (abs(metrics[[1]]$c1) >= abs(metrics[[1]]$c2)) j1 <- 1L else j2 <- 1L
  }
  aligned <- list(F1 = NULL, F2 = NULL); retention <- c(F1 = NA_real_, F2 = NA_real_)
  cross_abs <- c(F1 = NA_real_, F2 = NA_real_)
  for (d in 1:2) {
    j <- if (d == 1) j1 else j2
    if (!is.na(j)) {
      cc <- if (d == 1) metrics[[j]]$c1 else metrics[[j]]$c2
      aligned[[d]] <- L[, j] * ifelse(cc < 0, -1, 1)
      retention[d] <- abs(cc)
      cross_abs[d] <- abs(if (d == 1) metrics[[j]]$c2 else metrics[[j]]$c1)
    }
  }
  list(nf = nf, j = c(F1 = j1, F2 = j2), aligned = aligned, retention = retention,
       cross_abs = cross_abs, audit = if (length(factor_audit)) rbindlist(factor_audit) else data.frame())
}

matches <- list(Swedish = list(), HUNT = list()); audit_rows <- list()
for (cohort in c("Swedish", "HUNT")) {
  G <- if (cohort == "Swedish") G_SW else G_HU
  for (aid in all_ids) {
    z <- match_cell(loadings[[cohort]][[aid]], G, cohort, aid, pves[[cohort]][[aid]])
    matches[[cohort]][[aid]] <- z
    if (nrow(z$audit)) audit_rows[[length(audit_rows) + 1L]] <- z$audit
  }
}
factor_audit <- rbindlist(audit_rows, fill = TRUE)
factor_audit$one_to_one_selected <- FALSE
for (cohort in c("Swedish", "HUNT")) for (aid in all_ids) for (d in 1:2) {
  j <- matches[[cohort]][[aid]]$j[d]
  if (!is.na(j)) factor_audit$one_to_one_selected[factor_audit$annotation_id == aid & factor_audit$cohort == cohort &
    factor_audit$global_direction == paste0("F", d) & factor_audit$factor == paste0("F", j)] <- TRUE
}

subspace_metrics <- function(A, H) {
  if (!ncol(A) || !ncol(H)) return(c(full_cos1 = NA, full_cos2 = NA, full_angle1 = NA, full_angle2 = NA))
  qa <- qr.Q(qr(A)); qh <- qr.Q(qr(H)); sv <- pmin(1, svd(t(qa) %*% qh, nu = 0, nv = 0)$d)
  c(full_cos1 = sv[1], full_cos2 = if (length(sv) >= 2) sv[2] else NA,
    full_angle1 = acos(sv[1]) * 180 / pi, full_angle2 = if (length(sv) >= 2) acos(sv[2]) * 180 / pi else NA)
}

direction_rows <- list(); cell_rows <- list()
for (aid in all_ids) {
  swm <- matches$Swedish[[aid]]; hum <- matches$HUNT[[aid]]
  pref_sw <- if (all(is.finite(swm$retention))) (swm$retention[1] - swm$retention[2]) / sum(swm$retention) else NA_real_
  pref_hu <- if (all(is.finite(hum$retention))) (hum$retention[1] - hum$retention[2]) / sum(hum$retention) else NA_real_
  pref_ok <- if (is.finite(pref_sw) && is.finite(pref_hu)) (sign(pref_sw) == sign(pref_hu)) || (abs(pref_sw) < .1 && abs(pref_hu) < .1) else NA
  for (d in 1:2) {
    a <- swm$aligned[[d]]; h <- hum$aligned[[d]]; evaluable <- !is.null(a) && !is.null(h)
    if (evaluable) {
      oa <- rank_order(a, trait_order, species); oh <- rank_order(h, trait_order, species)
      rbo <- rbo_orders(oa, oh); o20 <- top_overlap(oa, oh, 20L); o40 <- top_overlap(oa, oh, 40L)
      dc20 <- direction_concordance(a, h, oa, oh, 20L)
      cs <- cosine(a, h); pe <- cor(a, h); sp <- cor(a, h, method = "spearman")
    } else { rbo <- o20 <- o40 <- dc20 <- cs <- pe <- sp <- NA_real_ }
    direction_rows[[length(direction_rows) + 1L]] <- data.frame(
      annotation_id = aid, analysis_family = if (aid %in% primary_ids) "PRIMARY_18" else "POST_STAGEA_EXPLORATORY",
      organ = ifelse(startsWith(aid, "gut__"), "gut", "brain"), direction = paste0("F", d),
      swedish_factor_count = swm$nf, hunt_factor_count = hum$nf,
      swedish_matched_factor = ifelse(is.na(swm$j[d]), NA, paste0("F", swm$j[d])),
      hunt_matched_factor = ifelse(is.na(hum$j[d]), NA, paste0("F", hum$j[d])), direction_evaluable = evaluable,
      swedish_global_retention_abs_cosine = swm$retention[d], hunt_global_retention_abs_cosine = hum$retention[d],
      swedish_global_specificity_margin = pmax(0, swm$retention[d] - swm$cross_abs[d]),
      hunt_global_specificity_margin = pmax(0, hum$retention[d] - hum$cross_abs[d]),
      crosscohort_signed_cosine = cs, crosscohort_pearson = pe, crosscohort_spearman = sp,
      crosscohort_RBO_p095 = rbo, crosscohort_top20_overlap = o20, crosscohort_top40_overlap = o40,
      common_top20_direction_concordance = dc20, swedish_preference = pref_sw, hunt_preference = pref_hu,
      preference_concordant = pref_ok, stringsAsFactors = FALSE)
  }
  sm <- subspace_metrics(loadings$Swedish[[aid]], loadings$HUNT[[aid]])
  ma <- if (all(!vapply(swm$aligned, is.null, logical(1)))) do.call(cbind, swm$aligned) else NULL
  mh <- if (all(!vapply(hum$aligned, is.null, logical(1)))) do.call(cbind, hum$aligned) else NULL
  ms <- if (!is.null(ma) && !is.null(mh)) subspace_metrics(ma, mh) else rep(NA_real_, 4)
  cell_rows[[length(cell_rows) + 1L]] <- data.frame(
    annotation_id = aid, analysis_family = if (aid %in% primary_ids) "PRIMARY_18" else "POST_STAGEA_EXPLORATORY",
    organ = ifelse(startsWith(aid, "gut__"), "gut", "brain"), swedish_factor_count = swm$nf, hunt_factor_count = hum$nf,
    factor_count_equal = swm$nf == hum$nf, swedish_both_directions = all(is.finite(swm$retention)),
    hunt_both_directions = all(is.finite(hum$retention)), both_directions_in_both_cohorts = all(is.finite(c(swm$retention, hum$retention))),
    swedish_preference = pref_sw, hunt_preference = pref_hu, preference_concordant = pref_ok,
    full_subspace_cosine_1 = sm[1], full_subspace_cosine_2 = sm[2], full_subspace_angle_1_degrees = sm[3], full_subspace_angle_2_degrees = sm[4],
    matched_global_subspace_cosine_1 = ms[1], matched_global_subspace_cosine_2 = ms[2],
    matched_global_subspace_angle_1_degrees = ms[3], matched_global_subspace_angle_2_degrees = ms[4], stringsAsFactors = FALSE)
}
directions <- rbindlist(direction_rows, fill = TRUE)
cells <- rbindlist(cell_rows, fill = TRUE)

# Reuse the exact Phase 7 primary information strata and generate a shared row-block permutation stream.
p7state <- readRDS(file.path(P7, "internal", "phase7_reproducibility_state.rds"))
strata <- p7state$primary_strata
stopifnot(length(strata) == 241L, identical(as.integer(table(strata)), c(26L, 23L, 31L, 29L, 26L, 25L, 25L, 31L, 25L)))
set.seed(SEED)
groups <- split(seq_len(241L), strata)
perm <- matrix(rep(seq_len(241L), each = B), nrow = B, ncol = 241L)
for (b in seq_len(B)) for (grp in groups) perm[b, grp] <- sample(grp, length(grp), replace = FALSE)
invperm <- matrix(NA_integer_, nrow = B, ncol = 241L)
for (b in seq_len(B)) invperm[b, perm[b, ]] <- seq_len(241L)

rbo_null <- function(swvec, huvec) {
  ord_a <- rank_order(swvec, trait_order, species); rank_a <- integer(241L); rank_a[ord_a] <- seq_len(241L)
  vals <- sort(unique(abs(huvec)), decreasing = TRUE)
  tie_groups <- lapply(vals, function(v) which(abs(huvec) == v))
  parts <- lapply(tie_groups, function(gr) {
    z <- invperm[, gr, drop = FALSE]
    if (length(gr) > 1L) t(apply(z, 1L, sort)) else z
  })
  new_order <- do.call(cbind, parts)
  rank_b <- matrix(NA_integer_, nrow = B, ncol = 241L)
  idx <- cbind(rep(seq_len(B), each = 241L), as.vector(t(new_order)))
  rank_b[idx] <- rep(seq_len(241L), times = B)
  m <- pmax(rank_b, matrix(rank_a, nrow = B, ncol = 241L, byrow = TRUE))
  w <- (1 - RBO_P) * RBO_P^(seq_len(241L) - 1L) / seq_len(241L)
  s <- rev(cumsum(rev(w)))
  out <- rowSums(matrix(s[m], nrow = B, ncol = 241L)) + RBO_P^241L
  # Verify the vectorized rank formula against direct row-block permutation.
  for (b in 1:3) {
    direct_h <- huvec[perm[b, ]]
    direct <- rbo_orders(ord_a, rank_order(direct_h, trait_order, species))
    stopifnot(abs(out[b] - direct) < 1e-12)
  }
  out
}

null_rows <- list()
for (i in seq_len(nrow(directions))) {
  row <- directions[i, ]; aid <- row$annotation_id; d <- as.integer(sub("F", "", row$direction))
  if (row$direction_evaluable) {
    nv <- rbo_null(matches$Swedish[[aid]]$aligned[[d]], matches$HUNT[[aid]]$aligned[[d]])
    obs <- row$crosscohort_RBO_p095
    ep <- (1 + sum(nv >= obs - 1e-15)) / (B + 1)
    null_rows[[i]] <- data.frame(annotation_id = aid, analysis_family = row$analysis_family, organ = row$organ,
      direction = row$direction, direction_evaluable = TRUE, permutations = B, seed = SEED,
      phase7_strata_sizes = paste(as.integer(table(strata)), collapse = ";"), observed_RBO_p095 = obs,
      null_mean = mean(nv), null_sd = sd(nv), null_median = median(nv), null_q95 = unname(quantile(nv, .95)),
      null_q99 = unname(quantile(nv, .99)), empirical_p = ep, stringsAsFactors = FALSE)
  } else {
    null_rows[[i]] <- data.frame(annotation_id = aid, analysis_family = row$analysis_family, organ = row$organ,
      direction = row$direction, direction_evaluable = FALSE, permutations = B, seed = SEED,
      phase7_strata_sizes = paste(as.integer(table(strata)), collapse = ";"), observed_RBO_p095 = NA,
      null_mean = NA, null_sd = NA, null_median = NA, null_q95 = NA, null_q99 = NA,
      empirical_p = 1, stringsAsFactors = FALSE)
  }
}
nulls <- rbindlist(null_rows)
nulls$BH_family <- ifelse(nulls$analysis_family == "PRIMARY_18", "PRIMARY_18x2_ALL_36", "EXPLORATORY_4x2_SEPARATE_8")
nulls$BH_q <- NA_real_
ix <- which(nulls$analysis_family == "PRIMARY_18"); stopifnot(length(ix) == 36L); nulls$BH_q[ix] <- p.adjust(nulls$empirical_p[ix], method = "BH")
ix2 <- which(nulls$analysis_family == "POST_STAGEA_EXPLORATORY"); stopifnot(length(ix2) == 8L); nulls$BH_q[ix2] <- p.adjust(nulls$empirical_p[ix2], method = "BH")
nulls$composition_validated_FDR_0_10 <- nulls$direction_evaluable & nulls$BH_q <= .10

directions <- merge(directions, nulls[, c("annotation_id", "direction", "empirical_p", "BH_family", "BH_q", "composition_validated_FDR_0_10")],
                    by = c("annotation_id", "direction"), all.x = TRUE, sort = FALSE)
directions$directionally_consistent <- directions$direction_evaluable & directions$crosscohort_signed_cosine > 0 &
  directions$common_top20_direction_concordance >= .50
directions$direction_supported <- directions$composition_validated_FDR_0_10 & directions$directionally_consistent

classify_cells <- function(family) {
  ids <- if (family == "PRIMARY_18") primary_ids else explore_ids
  out <- cells[match(ids, cells$annotation_id), ]
  out$F1_composition_validated <- out$F2_composition_validated <- FALSE
  out$F1_directionally_consistent <- out$F2_directionally_consistent <- FALSE
  out$F1_RBO <- out$F2_RBO <- out$F1_BH_q <- out$F2_BH_q <- NA_real_
  for (i in seq_len(nrow(out))) for (d in 1:2) {
    z <- directions[directions$annotation_id == out$annotation_id[i] & directions$direction == paste0("F", d), ]
    out[[paste0("F", d, "_composition_validated")]][i] <- z$composition_validated_FDR_0_10
    out[[paste0("F", d, "_directionally_consistent")]][i] <- z$directionally_consistent
    out[[paste0("F", d, "_RBO")]][i] <- z$crosscohort_RBO_p095
    out[[paste0("F", d, "_BH_q")]][i] <- z$BH_q
  }
  out$validation_status <- vapply(seq_len(nrow(out)), function(i) {
    if (!out$both_directions_in_both_cohorts[i] || is.na(out$preference_concordant[i])) return("NOT_EVALUABLE")
    supported_direction <- (out$F1_composition_validated[i] && out$F1_directionally_consistent[i]) ||
      (out$F2_composition_validated[i] && out$F2_directionally_consistent[i])
    if (supported_direction && out$preference_concordant[i]) "SUPPORTED" else "NOT_SUPPORTED"
  }, character(1))
  out
}
primary_results <- classify_cells("PRIMARY_18")
explore_results <- classify_cells("POST_STAGEA_EXPLORATORY")
explore_results$validation_status <- paste0("EXPLORATORY_", explore_results$validation_status)

wr(model_summary[model_summary$cohort == "HUNT", ], "STAGEB_HUNT_MODEL_SUMMARY.tsv")
wr(factor_audit, "STAGEB_HUNT_FACTOR_MATCHING_AUDIT.tsv")
wr(directions, "06_CELL_GFA_CROSSCOHORT_COMPOSITION.tsv")
wr(nulls[nulls$analysis_family == "PRIMARY_18", ], "STAGEB_HUNT_RBO_NULL_CALIBRATION.tsv")
wr(primary_results, "STAGEB_18_CELL_VALIDATION_RESULTS.tsv")
wr(explore_results, "STAGEB_POST_STAGEA_EXPLORATORY_RESULTS.tsv")

dir.create(file.path(P9, "internal"), showWarnings = FALSE)
saveRDS(list(seed = SEED, B = B, rbo_p = RBO_P, primary_strata = strata,
             primary_ids = primary_ids, exploratory_ids = explore_ids,
             model_summary = model_summary, directions = directions, nulls = nulls,
             primary_results = primary_results, exploratory_results = explore_results),
        file.path(P9, "internal", "stageb_validation_state.rds"), compress = "xz")

cat("Primary status:\n"); print(table(primary_results$validation_status, useNA = "ifany"))
cat("Exploratory status:\n"); print(table(explore_results$validation_status, useNA = "ifany"))
cat("Primary validated directions:", sum(directions$analysis_family == "PRIMARY_18" & directions$composition_validated_FDR_0_10), "/36\n")
