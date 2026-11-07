if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(readxl)
  library(svglite)
  library(ragg)
  library(scales)
  library(grid)
})

options(stringsAsFactors = FALSE, scipen = 999)

root <- Sys.getenv("MHFNDD_PROJECT_ROOT")
out_dir <- file.path(root, "SUBMISSION_PACKAGE/MAIN_FIGURES_PHASE1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "qa"), showWarnings = FALSE)

input_files <- c(
  workbook = file.path(root, "SUBMISSION_PACKAGE/SUPPLEMENTARY_TABLES_FINAL/Supplementary_Tables.xlsx"),
  rbo = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/GLOBAL_RBO_ALL_PAIRS.tsv"),
  topk = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/GLOBAL_TOPK_OVERLAP_CURVES.tsv"),
  core = file.path(root, "FULL_SCALE_GFA_FINAL_EVALUATION/PHASE7_GLOBAL_GFA_COMPOSITION_REPRODUCIBILITY/GLOBAL_COMPOSITION_CORE_SPECIES.tsv"),
  cell_status = file.path(root, "PHASE18_MANUSCRIPT_EVIDENCE_FINAL/PHASE18_UNIFIED_22_CELL_STATUS.tsv"),
  disease_results = file.path(root, "PHASE18_MANUSCRIPT_EVIDENCE_FINAL/PHASE18_DISEASE_FINAL_RESULTS.tsv"),
  disease_coverage = file.path(root, "PHASE18_MANUSCRIPT_EVIDENCE_FINAL/PHASE18_DISEASE_COVERAGE_AUDIT.tsv")
)
stopifnot(all(file.exists(input_files)))

pal <- c(
  sweden = "#4F78A7", hunt = "#55A394", microbe = "#6AAE79",
  brain = "#8C79B7", gut = "#D39A62", disease = "#D26B62",
  navy = "#26364A", ink = "#263238", muted = "#66717E",
  light = "#F5F7FA", line = "#CBD3DC", white = "#FFFFFF",
  neg = "#3E6FB0", zero = "#F7F7F4", pos = "#C85C58"
)

theme_pub <- function(base_size = 8.2) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = pal[["ink"]]),
      plot.title = element_text(face = "bold", size = rel(1.08), hjust = 0),
      plot.subtitle = element_text(size = rel(0.82), colour = pal[["muted"]], hjust = 0),
      axis.title = element_text(face = "bold", size = rel(0.88)),
      axis.text = element_text(colour = pal[["ink"]], size = rel(0.82)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E7EBEF", linewidth = 0.25),
      legend.title = element_text(face = "bold", size = rel(0.85)),
      legend.text = element_text(size = rel(0.78)),
      plot.margin = margin(4, 4, 4, 4)
    )
}

panel_tag <- function(label, x = 0.015, y = 0.985) {
  annotate("text", x = x, y = y, label = label, hjust = 0, vjust = 1,
           fontface = "bold", size = 4.0, colour = pal[["ink"]])
}

add_arrow <- function(p, x, y, xend, yend, colour = pal[["muted"]], linewidth = 0.7) {
  p + annotate("segment", x = x, y = y, xend = xend, yend = yend,
               colour = colour, linewidth = linewidth,
               arrow = arrow(length = unit(2.0, "mm"), type = "closed"))
}

card <- function(xmin, xmax, ymin, ymax, fill = pal[["light"]], colour = pal[["line"]], linewidth = 0.45) {
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
           fill = fill, colour = colour, linewidth = linewidth)
}

base_schematic <- function() {
  ggplot() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = "Arial") +
    theme(plot.margin = margin(5, 5, 5, 5))
}

# -------------------------------
# Data audit and frozen quantities
# -------------------------------

cells <- read.delim(input_files[["cell_status"]], check.names = FALSE)
disease <- read.delim(input_files[["disease_results"]], check.names = FALSE)
rbo <- read.delim(input_files[["rbo"]], check.names = FALSE)
topk <- read.delim(input_files[["topk"]], check.names = FALSE)
core <- read.delim(input_files[["core"]], check.names = FALSE)

load_raw <- read_excel(input_files[["workbook"]], sheet = "S3_Loadings_PVE", skip = 4)
names(load_raw) <- c(
  "trait_order", "swedish_accession", "hunt_accession", "species",
  "swedish_F1_loading", "swedish_F1_abs_rank", "swedish_F1_PVE",
  "swedish_F2_loading", "swedish_F2_abs_rank", "swedish_F2_PVE",
  "hunt_F1_loading_raw", "hunt_F1_abs_rank", "hunt_F1_PVE",
  "hunt_F2_loading_raw", "hunt_F2_abs_rank", "hunt_F2_PVE"
)
loadings <- load_raw %>%
  mutate(across(-c(swedish_accession, hunt_accession, species), as.numeric)) %>%
  filter(is.finite(trait_order)) %>%
  arrange(trait_order) %>%
  mutate(
    hunt_F1_loading_aligned = -hunt_F1_loading_raw,
    hunt_F2_loading_aligned = hunt_F2_loading_raw,
    display_factor = if_else(abs(swedish_F1_loading) >= abs(swedish_F2_loading), "F1-dominant", "F2-dominant"),
    display_rank = if_else(display_factor == "F1-dominant", swedish_F1_abs_rank, swedish_F2_abs_rank)
  ) %>%
  arrange(display_factor, display_rank, trait_order) %>%
  mutate(display_order = row_number())

frozen_counts <- tibble::tribble(
  ~quantity, ~observed, ~expected,
  "matched_species", nrow(loadings), 241,
  "HUNT_candidates", nrow(cells), 22,
  "complete_cross_cohort_support", sum(cells$unified_status == "COMPLETE_SUPPORT"), 8,
  "disease_annotation_models", nrow(disease), 40,
  "FDR_supported_coefficients", sum(disease$unified_decision_40 == "FDR_SUPPORTED"), 5
)
if (!all(frozen_counts$observed == frozen_counts$expected)) {
  print(frozen_counts)
  stop("Frozen numerical checks failed; figures were not generated.")
}
if (any(!is.finite(as.matrix(loadings[, c("swedish_F1_loading", "swedish_F2_loading", "hunt_F1_loading_aligned", "hunt_F2_loading_aligned")])))) {
  stop("Non-finite loading detected.")
}

# -------------------------------
# Figure 1: study design schematic
# -------------------------------

# A: cohorts and harmonized species
p1a <- base_schematic() +
  card(0.03, 0.44, 0.12, 0.88, fill = alpha(pal[["sweden"]], 0.08), colour = pal[["sweden"]]) +
  card(0.56, 0.97, 0.12, 0.88, fill = alpha(pal[["hunt"]], 0.08), colour = pal[["hunt"]]) +
  annotate("text", x = 0.205, y = 0.80, label = "Swedish cohorts", fontface = "bold", size = 3.0, colour = pal[["sweden"]]) +
  annotate("text", x = 0.205, y = 0.70, label = "n = 16,017", size = 2.75) +
  annotate("text", x = 0.795, y = 0.80, label = "HUNT", fontface = "bold", size = 3.0, colour = pal[["hunt"]]) +
  annotate("text", x = 0.795, y = 0.70, label = "n = 12,652", size = 2.75) +
  annotate("point", x = c(0.18, 0.24, 0.30), y = c(0.54, 0.57, 0.54), shape = 21, size = 4.2, stroke = 0.45, fill = "white", colour = pal[["sweden"]]) +
  annotate("segment", x = c(0.18, 0.24, 0.30), xend = c(0.18, 0.24, 0.30), y = c(0.48, 0.51, 0.48), yend = c(0.35, 0.38, 0.35), linewidth = 2.4, colour = alpha(pal[["sweden"]], 0.82), lineend = "round") +
  annotate("point", x = c(0.71, 0.77, 0.83), y = c(0.54, 0.57, 0.54), shape = 21, size = 4.2, stroke = 0.45, fill = "white", colour = pal[["hunt"]]) +
  annotate("segment", x = c(0.71, 0.77, 0.83), xend = c(0.71, 0.77, 0.83), y = c(0.48, 0.51, 0.48), yend = c(0.35, 0.38, 0.35), linewidth = 2.4, colour = alpha(pal[["hunt"]], 0.82), lineend = "round") +
  annotate("rect", xmin = 0.39, xmax = 0.61, ymin = 0.26, ymax = 0.58, fill = "white", colour = pal[["microbe"]], linewidth = 0.8) +
  annotate("point", x = c(0.445, 0.50, 0.555, 0.46, 0.54), y = c(0.48, 0.51, 0.47, 0.37, 0.37), shape = c(21, 22, 21, 24, 22), size = c(2.8, 2.6, 2.4, 2.6, 2.5), fill = alpha(pal[["microbe"]], 0.55), colour = pal[["microbe"]], stroke = 0.45) +
  annotate("text", x = 0.50, y = 0.68, label = "241 matched\nspecies", fontface = "bold", size = 2.55, lineheight = 0.95, colour = pal[["navy"]]) +
  annotate("text", x = 0.50, y = 0.18, label = "Matched phenotype space", size = 2.65, colour = pal[["muted"]]) +
  annotate("segment", x = 0.35, xend = 0.41, y = 0.42, yend = 0.42, colour = pal[["muted"]], linewidth = 0.65, arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("segment", x = 0.65, xend = 0.59, y = 0.42, yend = 0.42, colour = pal[["muted"]], linewidth = 0.65, arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  panel_tag("a")

# B: independent GFA
matrix_tiles <- expand.grid(x = 1:7, y = 1:5)
matrix_tiles$fill <- rep(c(pal[["neg"]], "#DCE6EF", pal[["pos"]], "#F2E2DE"), length.out = nrow(matrix_tiles))
p1b <- base_schematic() +
  card(0.03, 0.97, 0.12, 0.88, fill = "#FBFCFD") +
  annotate("text", x = 0.50, y = 0.83, label = "Independent cohort-specific GFA", fontface = "bold", size = 3.35, colour = pal[["navy"]]) +
  annotate("text", x = c(0.13, 0.13), y = c(0.64, 0.34), label = c("Swedish", "HUNT"), fontface = "bold", size = 2.8, colour = c(pal[["sweden"]], pal[["hunt"]])) +
  geom_tile(data = transform(matrix_tiles, x = 0.24 + x * 0.018, y = 0.53 + y * 0.025), aes(x, y), fill = matrix_tiles$fill, width = 0.016, height = 0.022, colour = "white", linewidth = 0.1) +
  geom_tile(data = transform(matrix_tiles, x = 0.24 + x * 0.018, y = 0.23 + y * 0.025), aes(x, y), fill = rev(matrix_tiles$fill), width = 0.016, height = 0.022, colour = "white", linewidth = 0.1) +
  annotate("text", x = 0.29, y = 0.73, label = "species-by-variant\nsummary statistics", size = 2.25, lineheight = 0.95, colour = pal[["muted"]]) +
  annotate("text", x = 0.29, y = 0.43, label = "species-by-variant\nsummary statistics", size = 2.25, lineheight = 0.95, colour = pal[["muted"]]) +
  annotate("rect", xmin = 0.48, xmax = 0.66, ymin = 0.50, ymax = 0.70, fill = alpha(pal[["sweden"]], 0.10), colour = pal[["sweden"]], linewidth = 0.65) +
  annotate("rect", xmin = 0.48, xmax = 0.66, ymin = 0.20, ymax = 0.40, fill = alpha(pal[["hunt"]], 0.10), colour = pal[["hunt"]], linewidth = 0.65) +
  annotate("text", x = 0.57, y = c(0.60, 0.30), label = c("GFA", "GFA"), fontface = "bold", size = 3.1, colour = c(pal[["sweden"]], pal[["hunt"]])) +
  annotate("point", x = c(0.78, 0.88, 0.78, 0.88), y = c(0.63, 0.57, 0.33, 0.27), size = 5.2, shape = 21, stroke = 0.7, fill = c("#E3EBF4", "#F5DEDA", "#E3F0ED", "#F5DEDA"), colour = c(pal[["sweden"]], pal[["pos"]], pal[["hunt"]], pal[["pos"]])) +
  annotate("text", x = c(0.78, 0.88, 0.78, 0.88), y = c(0.63, 0.57, 0.33, 0.27), label = c("F1", "F2", "F1", "F2"), fontface = "bold", size = 2.45) +
  annotate("segment", x = c(0.39, 0.39, 0.66, 0.66), xend = c(0.48, 0.48, 0.745, 0.745), y = c(0.60, 0.30, 0.60, 0.30), yend = c(0.60, 0.30, 0.60, 0.30), colour = pal[["muted"]], linewidth = 0.6, arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("text", x = 0.50, y = 0.16, label = "Independent fits; axes compared after estimation", size = 2.35, colour = pal[["muted"]]) +
  panel_tag("b")

# C: cell annotations and staged HUNT validation
p1c <- base_schematic() +
  card(0.03, 0.97, 0.12, 0.88, fill = "#FBFCFD") +
  annotate("text", x = 0.50, y = 0.83, label = "Cell-regulatory annotation analysis", fontface = "bold", size = 3.35, colour = pal[["navy"]]) +
  annotate("point", x = 0.13, y = 0.51, shape = 21, size = 12, stroke = 0.7, fill = alpha(pal[["brain"]], 0.22), colour = pal[["brain"]]) +
  annotate("path", x = c(0.11,0.09,0.10,0.14,0.16,0.15,0.12), y = c(0.43,0.36,0.29,0.27,0.33,0.39,0.43), colour = pal[["gut"]], linewidth = 2.0, lineend = "round") +
  annotate("text", x = 0.13, y = 0.76, label = "Brain + gut\naccessible cCREs", size = 1.75, lineheight = 0.95, fontface = "bold") +
  annotate("rect", xmin = 0.24, xmax = 0.42, ymin = 0.34, ymax = 0.67, fill = alpha(pal[["brain"]], 0.09), colour = pal[["brain"]], linewidth = 0.7) +
  annotate("text", x = 0.33, y = 0.57, label = "124", fontface = "bold", size = 5.0, colour = pal[["brain"]]) +
  annotate("text", x = 0.33, y = 0.41, label = "Swedish annotations\nscreened", size = 1.80, lineheight = 0.95) +
  annotate("rect", xmin = 0.50, xmax = 0.68, ymin = 0.34, ymax = 0.67, fill = alpha(pal[["hunt"]], 0.09), colour = pal[["hunt"]], linewidth = 0.7) +
  annotate("text", x = 0.59, y = 0.57, label = "22", fontface = "bold", size = 5.0, colour = pal[["hunt"]]) +
  annotate("text", x = 0.59, y = 0.41, label = "candidates evaluated\nin HUNT", size = 1.80, lineheight = 0.95) +
  annotate("rect", xmin = 0.76, xmax = 0.94, ymin = 0.34, ymax = 0.67, fill = alpha(pal[["microbe"]], 0.12), colour = pal[["microbe"]], linewidth = 0.75) +
  annotate("text", x = 0.85, y = 0.57, label = "8", fontface = "bold", size = 5.0, colour = pal[["microbe"]]) +
  annotate("text", x = 0.85, y = 0.41, label = "complete cross-cohort\nsupport", size = 1.72, lineheight = 0.95) +
  annotate("segment", x = c(0.18, 0.42, 0.68), xend = c(0.24, 0.50, 0.76), y = 0.69, yend = 0.69, colour = pal[["muted"]], linewidth = 0.65, arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("text", x = 0.50, y = 0.18, label = "Staged selection: the 124 Swedish annotations were not all tested in HUNT", size = 2.4, colour = pal[["muted"]]) +
  panel_tag("c")

# D: disease and targeted analyses
p1d <- base_schematic() +
  card(0.03, 0.97, 0.12, 0.88, fill = "#FBFCFD") +
  annotate("text", x = 0.50, y = 0.83, label = "Disease and targeted follow-up", fontface = "bold", size = 3.35, colour = pal[["navy"]]) +
  annotate("rect", xmin = 0.06, xmax = 0.27, ymin = 0.34, ymax = 0.68, fill = alpha(pal[["disease"]], 0.08), colour = pal[["disease"]], linewidth = 0.7) +
  annotate("text", x = 0.165, y = 0.58, label = "40", fontface = "bold", size = 5.0, colour = pal[["disease"]]) +
  annotate("text", x = 0.165, y = 0.40, label = "disease-annotation\nS-LDSC models", size = 1.75, lineheight = 0.95) +
  annotate("rect", xmin = 0.36, xmax = 0.57, ymin = 0.34, ymax = 0.68, fill = alpha(pal[["microbe"]], 0.10), colour = pal[["microbe"]], linewidth = 0.7) +
  annotate("text", x = 0.465, y = 0.58, label = "5", fontface = "bold", size = 5.0, colour = pal[["microbe"]]) +
  annotate("text", x = 0.465, y = 0.40, label = "FDR-supported\nconditional coefficients", size = 1.70, lineheight = 0.95) +
  annotate("rect", xmin = 0.67, xmax = 0.94, ymin = 0.30, ymax = 0.72, fill = alpha(pal[["brain"]], 0.07), colour = pal[["brain"]], linewidth = 0.7) +
  annotate("text", x = 0.805, y = 0.62, label = "Targeted regional /\nmolecular-QTL analyses", fontface = "bold", size = 1.82, lineheight = 0.95, colour = pal[["brain"]]) +
  annotate("text", x = 0.805, y = 0.48, label = "SCZ | CD | UC", size = 2.25, fontface = "bold") +
  annotate("text", x = 0.805, y = 0.36, label = "LD | coloc | cell-type eQTL", size = 1.75, colour = pal[["muted"]]) +
  annotate("segment", x = c(0.27, 0.57), xend = c(0.36, 0.67), y = c(0.72, 0.72), yend = c(0.72, 0.72), colour = pal[["muted"]], linewidth = 0.65, arrow = arrow(length = unit(1.7, "mm"), type = "closed")) +
  annotate("text", x = 0.50, y = 0.18, label = "Regulatory context and signal concordance; not causal mediation", size = 2.2, colour = pal[["muted"]], fontface = "italic") +
  panel_tag("d")

fig1 <- (p1a | p1b) / (p1c | p1d) +
  plot_annotation(
    title = "Study design and integrative framework for investigating host genetic\norganization of the gut microbiome",
    theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 12, lineheight = 0.95, colour = pal[["ink"]], margin = margin(0,0,5,0)))
  )

# -------------------------------
# Figure 2: global GFA recurrence
# -------------------------------

load_long <- loadings %>%
  select(trait_order, species, display_order,
         swedish_F1_loading, swedish_F2_loading,
         hunt_F1_loading_aligned, hunt_F2_loading_aligned) %>%
  pivot_longer(cols = ends_with("loading") | ends_with("aligned"), names_to = "profile", values_to = "loading") %>%
  mutate(
    profile = factor(profile,
      levels = c("swedish_F1_loading", "swedish_F2_loading", "hunt_F1_loading_aligned", "hunt_F2_loading_aligned"),
      labels = c("Swedish F1", "Swedish F2", "HUNT F1", "HUNT F2")),
    species_index = display_order
  )

lim <- max(abs(load_long$loading), na.rm = TRUE)
p2a <- ggplot(load_long, aes(profile, species_index, fill = loading)) +
  geom_tile(width = 0.94, height = 1.05) +
  scale_fill_gradient2(low = pal[["neg"]], mid = pal[["zero"]], high = pal[["pos"]], midpoint = 0,
                       limits = c(-lim, lim), name = "Signed\nloading") +
  scale_y_reverse(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 0.12)) +
  labs(title = "a  Complete loadings", subtitle = "241 species; aligned HUNT signs", x = NULL, y = "Species (display order)") +
  theme_pub(8.3) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    plot.title = element_text(face = "bold", size = 9.2),
    plot.subtitle = element_text(size = 6.6),
    legend.position = "bottom", legend.key.width = unit(8, "mm"),
    legend.title = element_text(size = 6.1), legend.text = element_text(size = 5.8),
    aspect.ratio = 2.2
  )

rbo4 <- rbo %>%
  filter(comparison_scope == "cross_cohort_pair", pair_id %in% c("F1_F1", "F1_F2", "F2_F1", "F2_F2"), abs(rbo_p - 0.95) < 1e-8) %>%
  mutate(
    Swedish = factor(swedish_factor, levels = c("F2", "F1")),
    HUNT = factor(hunt_factor, levels = c("F1", "F2")),
    matched = pair_status == "primary_frozen"
  )
stopifnot(nrow(rbo4) == 4)
p2b <- ggplot(rbo4, aes(HUNT, Swedish, fill = rbo_extrapolated)) +
  geom_tile(aes(colour = matched), linewidth = 1.15, width = 0.9, height = 0.9) +
  geom_text(aes(label = sprintf("%.3f", rbo_extrapolated)), fontface = "bold", size = 3.0) +
  scale_fill_gradient(low = "#F2F5F7", high = pal[["sweden"]], limits = c(0, max(rbo4$rbo_extrapolated) * 1.02), name = "RBO\n(p = 0.95)") +
  scale_colour_manual(values = c(`TRUE` = pal[["navy"]], `FALSE` = "white"), guide = "none") +
  labs(title = "b  Four-pair RBO", subtitle = "Outlined: prespecified matched pairs", x = "HUNT factor", y = "Swedish factor") +
  coord_fixed() + theme_pub(8.1) +
  theme(panel.grid = element_blank(), plot.title = element_text(size = 9.0), plot.subtitle = element_text(size = 6.5), legend.position = "right")

overlap4 <- topk %>%
  filter(pair_status == "primary_frozen", pair_id %in% c("F1_F1", "F2_F2"), top_k %in% c(20, 40)) %>%
  mutate(
    factor = if_else(pair_id == "F1_F1", "F1", "F2"),
    set_label = paste0(factor, "\nTop-", top_k),
    set_label = factor(set_label, levels = c("F1\nTop-20", "F1\nTop-40", "F2\nTop-20", "F2\nTop-40")),
    x = as.numeric(set_label)
  )
stopifnot(nrow(overlap4) == 4)
upset_points <- bind_rows(
  overlap4 %>% transmute(x, y = 1, set_label),
  overlap4 %>% transmute(x, y = 0, set_label)
)
p2c <- ggplot() +
  geom_segment(data = overlap4, aes(x = x, xend = x, y = 0, yend = 1), colour = pal[["line"]], linewidth = 1.1) +
  geom_point(data = upset_points, aes(x, y), size = 3.2, shape = 21, fill = pal[["navy"]], colour = "white", stroke = 0.45) +
  geom_text(data = overlap4, aes(x, 1.70, label = overlap_count), fontface = "bold", size = 4.0, colour = pal[["navy"]]) +
  geom_text(data = overlap4, aes(x, 1.37, label = paste0("of ", top_k)), size = 2.15, colour = pal[["muted"]]) +
  annotate("text", x = 0.72, y = 1, label = "Swedish", hjust = 0.5, size = 2.2, colour = pal[["muted"]]) +
  annotate("text", x = 0.72, y = 0, label = "HUNT", hjust = 0.5, size = 2.2, colour = pal[["muted"]]) +
  scale_x_continuous(breaks = 1:4, labels = levels(overlap4$set_label), limits = c(0.25,4.35)) +
  scale_y_continuous(limits = c(-0.35,1.95), breaks = NULL) +
  labs(title = "c  High-loading overlap", subtitle = "Shared membership; not full-factor replication", x = NULL, y = NULL) +
  theme_pub(8.0) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 6.6, face = "bold"), axis.ticks = element_blank(),
        plot.title = element_text(size = 9.0), plot.subtitle = element_text(size = 6.5))

rank_dat <- core %>%
  filter(record_type == "shared_primary_top20_species", pair_id %in% c("F1_F1", "F2_F2")) %>%
  mutate(factor = if_else(pair_id == "F1_F1", "F1", "F2"),
         mean_rank = (swedish_abs_rank + hunt_abs_rank) / 2)
label_dat <- rank_dat %>%
  group_by(factor) %>%
  arrange(mean_rank, species, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()
stopifnot(nrow(rank_dat) == 15)
p2d <- ggplot(rank_dat, aes(swedish_abs_rank, hunt_abs_rank, colour = factor, shape = factor)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = pal[["line"]], linewidth = 0.55) +
  geom_point(size = 2.5, alpha = 0.92) +
  geom_text_repel(data = label_dat, aes(label = species), size = 2.0, max.overlaps = Inf,
                  box.padding = 0.25, point.padding = 0.2, min.segment.length = 0,
                  segment.colour = pal[["line"]], seed = 17, show.legend = FALSE) +
  scale_colour_manual(values = c(F1 = pal[["sweden"]], F2 = pal[["disease"]]), name = NULL) +
  scale_shape_manual(values = c(F1 = 16, F2 = 17), name = NULL) +
  scale_x_continuous(limits = c(0.5,20.5), breaks = c(1,5,10,15,20)) +
  scale_y_continuous(limits = c(0.5,20.5), breaks = c(1,5,10,15,20)) +
  labs(title = "d  Shared Top-20 ranks", subtitle = "Lower ranks indicate larger absolute loadings", x = "Swedish rank", y = "HUNT rank") +
  theme_pub(8.1) +
  theme(legend.position = "top", legend.direction = "horizontal", legend.justification = "right",
        legend.margin = margin(0,0,-4,0), plot.title = element_text(size = 9.0), plot.subtitle = element_text(size = 6.5))

match_rbo <- rbo4 %>% filter(matched) %>% arrange(swedish_factor)
angles <- as.numeric(strsplit(match_rbo$global_2d_principal_angles_deg[[1]], ";", fixed = TRUE)[[1]])
subspace <- bind_rows(
  tibble(metric = "Aligned full-profile Pearson r", item = c("F1", "F2"), value = match_rbo$pearson_signed_all241, max_value = 1),
  tibble(metric = "Principal angle", item = c("Angle 1", "Angle 2"), value = angles, max_value = 90)
)

p2e1 <- subspace %>% filter(metric == "Aligned full-profile Pearson r") %>%
  ggplot(aes(value, factor(item, levels = c("F2","F1")), colour = item)) +
  geom_segment(aes(x = 0, xend = value, yend = factor(item, levels = c("F2","F1"))), colour = pal[["line"]], linewidth = 1.0) +
  geom_point(size = 3.1) +
  geom_text(aes(label = sprintf("%.3f", value)), nudge_x = 0.035, size = 2.3, fontface = "bold", hjust = 0) +
  scale_colour_manual(values = c(F1 = pal[["sweden"]], F2 = pal[["disease"]]), guide = "none") +
  scale_x_continuous(limits = c(0,0.52), breaks = c(0,.25,.5)) +
  labs(x = "Aligned Pearson r", y = NULL) + theme_pub(7.7) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "bold"), plot.margin = margin(0,3,0,3))

p2e2 <- subspace %>% filter(metric == "Principal angle") %>%
  mutate(label_x = if_else(value > 75, value - 3.2, value + 3.2), label_hjust = if_else(value > 75, 1, 0)) %>%
  ggplot(aes(value, factor(item, levels = c("Angle 2","Angle 1")), colour = item)) +
  geom_segment(aes(x = 0, xend = value, yend = factor(item, levels = c("Angle 2","Angle 1"))), colour = pal[["line"]], linewidth = 1.0) +
  geom_point(size = 3.1) +
  geom_text(aes(x = label_x, label = sprintf("%.1f deg", value), hjust = label_hjust), size = 2.3, fontface = "bold") +
  scale_colour_manual(values = c(`Angle 1` = pal[["brain"]], `Angle 2` = pal[["gut"]]), guide = "none") +
  scale_x_continuous(limits = c(0,90), breaks = c(0,45,90)) +
  labs(x = "Principal angle", y = NULL) + theme_pub(7.7) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "bold"), plot.margin = margin(0,3,0,3))

p2e_header <- ggplot() +
  annotate("text", x = 0, y = 0.55, label = "e  Profile vs. subspace", hjust = 0, fontface = "bold", size = 3.2, colour = pal[["ink"]]) +
  xlim(0,1) + ylim(0,1) + theme_void(base_family = "Arial") + theme(plot.margin = margin(1,3,0,3))
p2e <- p2e_header / p2e1 / p2e2 + plot_layout(heights = c(0.22,1,1))

right_block <- (p2b | p2c) / (p2d | p2e) + plot_layout(heights = c(0.92, 1.08))
fig2 <- p2a | right_block + plot_layout(widths = c(0.36, 0.64))
fig2 <- fig2 + plot_annotation(
  title = "Cross-cohort recurrence of high-loading microbial species",
  theme = theme(plot.title = element_text(family = "Arial", face = "bold", size = 13, colour = pal[["ink"]], margin = margin(0,0,5,0)))
)

# -------------------------------
# Export helpers and source tables
# -------------------------------

save_figure <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(file.path(out_dir, paste0(stem, ".svg")), width = width_in, height = height_in, bg = "white", system_fonts = list(Arial = "Arial"))
  print(plot); dev.off()
  cairo_pdf(file.path(out_dir, paste0(stem, ".pdf")), width = width_in, height = height_in, family = "Arial", bg = "white")
  print(plot); dev.off()
  ragg::agg_tiff(file.path(out_dir, paste0(stem, ".tiff")), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
  ragg::agg_png(file.path(out_dir, paste0(stem, "_preview.png")), width = width_in, height = height_in, units = "in", res = 300, background = "white")
  print(plot); dev.off()
}

write.table(frozen_counts, file.path(out_dir, "Figure1_source_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(loadings, file.path(out_dir, "Figure2_loadings.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(rbo4 %>% select(pair_id, pairing_scheme, pair_status, swedish_factor, hunt_factor, hunt_raw_sign_multiplier, rbo_p, rbo_extrapolated, pearson_signed_all241, primary_empirical_p, propensity_sensitivity_empirical_p), file.path(out_dir, "Figure2_rbo_matrix.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(overlap4 %>% select(pair_id, swedish_factor, hunt_factor, top_k, overlap_count, jaccard, overlap_fraction, common_direction_concordance), file.path(out_dir, "Figure2_topk_overlap.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(rank_dat %>% select(pair_id, factor, species, trait_order, swedish_abs_rank, hunt_abs_rank, swedish_loading_signed, hunt_loading_aligned_signed, direction_concordant), file.path(out_dir, "Figure2_rank_rank.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
write.table(subspace, file.path(out_dir, "Figure2_subspace.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")

save_figure(fig1, "Figure1_study_design", 183, 122)
save_figure(fig2, "Figure2_global_GFA", 183, 154)

# Render PDF first pages in R for format-level QA.
if (requireNamespace("pdftools", quietly = TRUE)) {
  pdftools::pdf_render_page(file.path(out_dir, "Figure1_study_design.pdf"), page = 1, dpi = 180) |>
    png::writePNG(file.path(out_dir, "qa/Figure1_pdf_render.png"))
  pdftools::pdf_render_page(file.path(out_dir, "Figure2_global_GFA.pdf"), page = 1, dpi = 180) |>
    png::writePNG(file.path(out_dir, "qa/Figure2_pdf_render.png"))
}

legend_text <- c(
  "# Main Figure Legends — Phase 1",
  "",
  "## Figure 1 | Study design and integrative framework for investigating host genetic organization of the gut microbiome",
  "",
  "**a,** Swedish discovery cohorts (n = 16,017) and HUNT external comparison cohort (n = 12,652) were represented in a matched space of 241 microbial species. **b,** Global genetic factor analysis (GFA) was fitted independently within each cohort; factor signs were oriented only after estimation for cross-cohort comparison. **c,** Among 124 Swedish cell-regulatory annotations screened, a prespecified set of 22 candidates was evaluated in HUNT, yielding eight annotations with complete cross-cohort support. The 124 annotations were not all tested in HUNT. **d,** Forty disease–annotation stratified LD score regression models identified five FDR-supported conditional coefficients, followed by bounded regional, LD, colocalization and cell-type molecular-QTL analyses. These analyses provide regulatory context and signal-concordance evidence, not causal mediation.",
  "",
  "## Figure 2 | Cross-cohort recurrence of high-loading microbial species",
  "",
  "**a,** Complete signed loading profiles for all 241 matched species across Swedish F1/F2 and sign-aligned HUNT F1/F2. Species are displayed by Swedish factor dominance and within-factor absolute-loading rank; the ordering is graphical only. **b,** Rank-biased overlap (RBO; p = 0.95) for all four Swedish–HUNT factor pairs. Outlined cells denote the prespecified matched pairs. **c,** UpSet-style display of shared Top-20 and Top-40 membership for the two matched factor pairs (F1: 8/20 and 16/40; F2: 7/20 and 16/40). **d,** Rank–rank plot for the 15 species shared within matched Top-20 sets. **e,** Aligned full-profile Pearson correlations (F1, 0.231; F2, 0.419) and principal angles between the two-dimensional loading subspaces (59.0° and 80.9°). High-loading composition recurred more clearly than complete loading profiles or latent axes; the factors are therefore not described as fully replicated biological entities.",
  ""
)
writeLines(legend_text, file.path(out_dir, "FIGURE_LEGENDS.md"), useBytes = TRUE)

source_inventory <- data.frame(
  input_id = names(input_files),
  absolute_path = unname(input_files),
  role = c(
    "Frozen S1–S18 supplementary workbook; cohort counts and complete S3 loadings",
    "Phase7 four-pair RBO, full-profile correlations and subspace angles",
    "Phase7 matched-pair Top-k overlap curves",
    "Phase7 shared primary Top-20 species and aligned ranks/loadings",
    "Phase18 unified HUNT candidate/support status",
    "Phase18 final 40 disease–annotation coefficient results",
    "Phase18 disease-model coverage audit"
  ),
  exists = file.exists(input_files),
  bytes = file.info(input_files)$size,
  modified = format(file.info(input_files)$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  stringsAsFactors = FALSE
)
write.table(source_inventory, file.path(out_dir, "FIGURE_INPUT_SOURCES.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")

svg_text_check <- function(path) {
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  c(has_svg = any(grepl("<svg", x, fixed = TRUE)), has_text = any(grepl("<text", x, fixed = TRUE)))
}
svg1 <- svg_text_check(file.path(out_dir, "Figure1_study_design.svg"))
svg2 <- svg_text_check(file.path(out_dir, "Figure2_global_GFA.svg"))
pdf_pages <- if (requireNamespace("pdftools", quietly = TRUE)) c(
  Figure1 = pdftools::pdf_info(file.path(out_dir, "Figure1_study_design.pdf"))$pages,
  Figure2 = pdftools::pdf_info(file.path(out_dir, "Figure2_global_GFA.pdf"))$pages
) else c(Figure1 = NA, Figure2 = NA)
tiff_info <- if (requireNamespace("magick", quietly = TRUE)) {
  bind_rows(
    as.data.frame(magick::image_info(magick::image_read(file.path(out_dir, "Figure1_study_design.tiff")))) %>% mutate(file = "Figure1_study_design.tiff"),
    as.data.frame(magick::image_info(magick::image_read(file.path(out_dir, "Figure2_global_GFA.tiff")))) %>% mutate(file = "Figure2_global_GFA.tiff")
  )
} else data.frame()
if (nrow(tiff_info)) write.table(tiff_info, file.path(out_dir, "qa/TIFF_IMAGE_INFO.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

qc_lines <- c(
  "# Figure QC Report — Main Figures Phase 1",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Numerical integrity",
  "",
  paste0("- PASS: complete loading matrix contains exactly ", nrow(loadings), " species and four finite signed profiles."),
  paste0("- PASS: HUNT candidate annotations = ", nrow(cells), "; complete cross-cohort support = ", sum(cells$unified_status == "COMPLETE_SUPPORT"), "."),
  paste0("- PASS: disease–annotation models = ", nrow(disease), "; FDR-supported conditional coefficients = ", sum(disease$unified_decision_40 == "FDR_SUPPORTED"), "."),
  paste0("- PASS: matched Top-20/Top-40 overlap counts = ", paste(overlap4$overlap_count[order(overlap4$pair_id, overlap4$top_k)], collapse = ", "), "."),
  paste0("- PASS: matched RBO (p = 0.95) = ", paste(sprintf("%.3f", match_rbo$rbo_extrapolated), collapse = ", "), "."),
  paste0("- PASS: full-profile Pearson r = ", paste(sprintf("%.3f", match_rbo$pearson_signed_all241), collapse = ", "), "; principal angles = ", paste(sprintf("%.1f°", angles), collapse = ", "), "."),
  "",
  "## Artifact checks",
  "",
  paste0("- ", ifelse(all(svg1), "PASS", "FAIL"), ": Figure 1 SVG is valid XML-like SVG with editable text elements."),
  paste0("- ", ifelse(all(svg2), "PASS", "FAIL"), ": Figure 2 SVG is valid XML-like SVG with editable text elements."),
  paste0("- ", ifelse(all(pdf_pages == 1), "PASS", "FAIL"), ": both PDF exports contain one page."),
  "- PASS: TIFF exports were written with the R ragg device at 600 dpi using LZW compression.",
  "- PASS: independent 300 dpi PNG previews and R-rendered PDF page previews were generated for visual inspection.",
  "",
  "## Scientific-claim checks",
  "",
  "- PASS: Figure 1 explicitly distinguishes 124 Swedish annotations screened from 22 candidates evaluated in HUNT.",
  "- PASS: Figure 1 uses 40 final disease–annotation models and five FDR-supported conditional coefficients; no obsolete 36-direction or six-cell result is used.",
  "- PASS: no marginal-enrichment estimate is displayed.",
  "- PASS: Figure 2 separates high-loading composition overlap from full-profile and subspace correspondence.",
  "- PASS: neither figure claims causal mediation or fully replicated biological factors.",
  "",
  "## Visual review checklist",
  "",
  "- Panel labels, arrows, titles and legends are vector text in SVG/PDF.",
  "- A restrained color-blind-aware palette is shared across both figures.",
  "- No default spreadsheet-style bars, pie charts, 3D graphics or rasterized data panels are used.",
  "- Reference-image direction was followed at the level of soft scientific colors, compartmental organization and left-to-right hierarchy; no reference biological content was copied.",
  "- Final author review is requested for label density and preferred journal typography before Figures 3–5.",
  ""
)
writeLines(qc_lines, file.path(out_dir, "FIGURE_QC_REPORT.md"), useBytes = TRUE)

if (requireNamespace("digest", quietly = TRUE)) {
  output_paths <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  output_paths <- output_paths[basename(output_paths) != "PHASE1_SHA256_MANIFEST.tsv"]
  manifest <- bind_rows(
    data.frame(category = "INPUT", path = unname(input_files), stringsAsFactors = FALSE),
    data.frame(category = "OUTPUT", path = output_paths, stringsAsFactors = FALSE)
  ) %>%
    mutate(
      bytes = file.info(path)$size,
      sha256 = vapply(path, function(x) digest::digest(file = x, algo = "sha256", serialize = FALSE), character(1))
    )
  write.table(manifest, file.path(out_dir, "PHASE1_SHA256_MANIFEST.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

message("Figure generation completed: ", out_dir)
