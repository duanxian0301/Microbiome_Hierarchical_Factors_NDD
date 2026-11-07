# Analysis and output map

| Module | Principal public scripts | Public derived outputs | Manuscript use |
|---|---|---|---|
| GWAS harmonization and common universe | `scripts/preprocessing/01_harmonize_trait.py`, `05_ld_prune.py`, `08_build_common_training_indices.py`, `10_build_r241_ldsc_inputs.py`, `11_run_official_r241_ldsc.R` | Parameters and accessions in `metadata/` | Methods; inputs to global and cell GFA |
| Global cohort-specific GFA | `scripts/global_gfa/15_fit_swedish241_rldsc.R`, `21_fit_hunt241_rldsc.R` | Figure 2 loading data | Global factor results |
| Cross-cohort comparison | `scripts/cross_cohort_validation/25_extract_and_compare_primary_factors.R` | Figure 2 RBO, overlap, rank, and correspondence data | Cross-cohort recurrence and subspace correspondence |
| Cell-regulatory validation | `scripts/cross_cohort_validation/21_stageb_crosscohort_validation.R` | `PHASE18_HUNT_44_TEST_FINAL.tsv`, `PHASE18_UNIFIED_22_CELL_STATUS.tsv`, Figure 4 source data | Cell-regulatory context |
| Disease S-LDSC | `scripts/disease_association/run_primary_phase12a.py`, `parse_primary_phase12a.py`, `run_phase17.py`, `finalize_phase17.py` | `PHASE18_DISEASE_FINAL_RESULTS.tsv`, `PHASE18_DISEASE_UNIFIED_BH.tsv`, Figure 4 source data | Disease-associated annotation analysis |
| Targeted regional evidence | `scripts/regional_analysis/run_phase13b.py`, `run_phase14b.py`, `run_coloc_abf.R` | `PHASE18_LOCUS_MOLECULAR_EVIDENCE.tsv`, Figure 5 source data | SCZ, CD, and UC targeted loci |
| SCZ chr3 concordance audit | `scripts/regional_analysis/run_phase16_pairwise_coloc.R` | SCZ rows in final locus evidence and Figure 5 source data | Pairwise coloc sensitivity and direction audit |
| Final evidence integration | `scripts/disease_association/build_phase18.py` | Final Phase18 tables | Results evidence harmonization |
| Figure generation | `scripts/figures/` | `results/figure_source_data/` | Figures 1–5 |

The main and supplementary manuscript documents are not included. The supplementary-table source crosswalk and index are retained in `results/supplementary/`.
