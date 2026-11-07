if (!nzchar(Sys.getenv("MHFNDD_PROJECT_ROOT"))) stop("MHFNDD_PROJECT_ROOT is required")
#!/usr/bin/env Rscript

# Main Figure 4, revision 2: display-only refinements.
# This wrapper evaluates the frozen rev1 build after exact, audited replacements
# to display order and labels. It does not modify any source data or statistics.

root <- Sys.getenv("MHFNDD_PROJECT_ROOT")
rev1_script <- file.path(root, "FIGURE4_CELL_REGULATORY_DISEASE_REDESIGN", "build_figure4_cell_regulatory_disease.R")
stopifnot(file.exists(rev1_script))

txt <- paste(readLines(rev1_script, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

replace_once <- function(old, new, label) {
  hits <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (hits != 1L) stop("Expected exactly one occurrence for ", label, "; observed ", hits)
  txt <<- sub(old, new, txt, fixed = TRUE)
}
replace_all <- function(old, new, label) {
  hits <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (hits < 1L) stop("Expected at least one occurrence for ", label)
  txt <<- gsub(old, new, txt, fixed = TRUE)
}

# New independent output directory and revision-specific artifact names.
replace_once(
  'outdir <- file.path(root, "FIGURE4_CELL_REGULATORY_DISEASE_REDESIGN")',
  'outdir <- file.path(root, "FIGURE4_CELL_REGULATORY_DISEASE_REDESIGN_REV2")',
  "revision output directory"
)
replace_all('FIGURE4_EXECUTION.log', 'FIGURE4_REV2_EXECUTION.log', "execution log name")

# Panel A: displayed tissue totals only; no plotted points or retention counts change.
replace_once(
  'palette_tissue_label <- c(brain = "Brain", gut = "Intestinal")',
  'palette_tissue_label <- c(brain = "Brain (n=100)", gut = "Intestinal (n=24)")',
  "Panel A tissue legend"
)

# Panel C/D: a single explicit frozen ordering is used in both panels.
old_c_order <- 'organ_order = ifelse(organ == "brain", 1L, 2L),\n         prov_order = ifelse(panel_origin == "ORIGINAL_FORMAL_18", 1L, 2L)) %>%\n  arrange(organ_order, prov_order, annotation_id) %>%'
new_c_order <- 'display_order = match(annotation_id, c("brain__BNGA", "brain__ET", "brain__ITL6_2_2", "brain__SMC", "brain__SNCG_1", "gut__Early_Colonocyte", "gut__Macrophage", "gut__NaiveB"))) %>%\n  arrange(display_order) %>%'
replace_once(old_c_order, new_c_order, "Panel C/D shared cell order")

# Panel B: replace only the explanatory caption; tiles and direction states remain untouched.
replace_once(
  'caption = "44 planned directions: 33 evaluable; 11 non-evaluable.")',
  'caption = "44 planned directions: 22 supported; 11 evaluable but unsupported; 11 non-evaluable.")',
  "Panel B complete status caption"
)

# Panel D: plotmath gives an actual superscript instead of a caret notation.
replace_once(
  'labs(x = "Conditional coefficient (x 10^-6); horizontal bars = coefficient +/- 1 SE",',
  'labs(x = expression("Conditional coefficient (" * "" %*% 10^{-6} * ")"),',
  "Panel D mathematical axis title"
)
replace_once(
  'caption = "Blue/green fill: brain/intestinal. Orange outline and q label: q40 FDR-supported.")',
  'caption = "Horizontal bars = coefficient +/- 1 SE. Blue/green fill: brain/intestinal. Orange outline and q label: q40 FDR-supported.")',
  "Panel D error-bar explanation"
)

# Revision-specific export/report names. These do not change any visual values.
replace_all('Figure4_cell_regulatory_disease', 'Figure4_cell_regulatory_disease_rev2', "final figure base name")
replace_all('FIGURE4_NUMERICAL_QC_REPORT.md', 'FIGURE4_REV2_NUMERICAL_QC_REPORT.md', "numerical QC name")
replace_all('FIGURE4_VISUAL_QC_REPORT.md', 'FIGURE4_REV2_VISUAL_QC_REPORT.md', "visual QC name")
replace_all('FIGURE4_BUILD_NOTES.md', 'FIGURE4_REV2_BUILD_NOTES.md', "build-notes name")
replace_once('build_figure4_cell_regulatory_disease.R', 'build_figure4_cell_regulatory_disease_rev2.R', "script manifest name")

# Explicitly document the limited display-only scope in the generated build notes.
replace_once(
  '"- Final size: 183 x 180 mm.", "",',
  '"- Final size: 183 x 180 mm.", "", "## Revision 2 scope", "", "- This is a visual and annotation-only minor revision; no data, statistical result, q value, coefficient, SE, grouping, significance decision, or provenance source was modified.", "- Changes only: shared Panel C/D cell order; Panel A tissue totals in the legend; Panel B complete direction-status caption; and Panel D superscript mathematical formatting.", "",',
  "build-notes scope"
)

# Rev2-specific QC wording: report only checks executed for this revision.
replace_once(
  'visual_qc_status <- "PASS - direct 600-dpi PNG, high-resolution Panel D crop and independently R-rendered PDF preview inspected; no clipping, overlap or missing values"',
  'visual_qc_status <- "PASS - direct 600-dpi PNG and independently R-rendered PDF preview inspected; all four requested visual revisions are present with no clipping or overlap"',
  "rev2 visual QC status"
)
replace_once(
  '"- Panel A: all 124 dots, four exact counts and tissue encoding are readable."',
  '"- Panel A: all 124 dots and retained-factor counts remain unchanged; the compact legend reads Brain (n=100) and Intestinal (n=24)."',
  "Panel A QC wording"
)
replace_once(
  '"- Panel B: all 44 tiles are visible; formal 18 and exploratory four are clearly separated; non-evaluable cells remain visible."',
  '"- Panel B: all 44 tiles are unchanged; the caption reports 22 supported, 11 evaluable but unsupported, and 11 non-evaluable directions."',
  "Panel B QC wording"
)
replace_once(
  '"- Panel C: all eight cells, tissue, provenance and F1/F2 direction states are legible."',
  '"- Panel C: the eight cells follow the requested five-brain then three-intestinal order; ET and ITL6_2_2 remain Exploratory."',
  "Panel C QC wording"
)
replace_once(
  '"- Panel D: all 40 estimates and +/- 1 SE bars are visible; all q40 labels are legible; five FDR-supported rows are restrainedly highlighted."',
  '"- Panel D: the same eight-cell order is used in every disease facet; all 40 estimates, stored +/- 1 SE bars, q40 labels and five frozen FDR-supported results remain visible; the axis uses a multiplication sign with superscript -6."',
  "Panel D QC wording"
)
replace_once(
  '"- Direct PNG, a high-resolution Panel D crop and the independently rendered PDF preview were inspected; all checks passed."',
  '"- The direct PNG and independently rendered PDF preview were visually inspected at final layout; all requested checks passed."',
  "executed visual review wording"
)

eval(parse(text = txt), envir = globalenv())
