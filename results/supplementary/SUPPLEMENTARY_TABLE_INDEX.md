# Supplementary Table Index

This submission package contains **18 numbered supplementary tables**. The main workbook has **38 data worksheets** plus the Guide worksheet. The complete S14B regional screening dataset is supplied separately as **Supplementary Data 1**. Lettered worksheets and the separate S14B attachment remain components of the same numbered table.

| Table | Worksheets | Total data rows | Scientific scope |
|---|---|---:|---|
| S1 | `S1A_Cohorts`, `S1B_Public_Resources` | 16 | Microbiome cohorts and variant resources; Public datasets, phenotypes and accessions |
| S2 | `S2_Species_Heritability` | 241 | Species identity and SNP heritability |
| S3 | `S3_Loadings_PVE` | 241 | Complete factor loadings, ranks and species-level PVE |
| S4 | `S4A_Factor_Correlations`, `S4B_Factor_Matching` | 18 | Cross-cohort factor correlations; Optimal factor matching and sign alignment |
| S5 | `S5A_TopK_Overlap`, `S5B_Shared_Species` | 51 | Top-20 and Top-40 species overlap; Species shared across matched high-loading sets |
| S6 | `S6A_RBO_Permutation`, `S6B_TopK_Curves`, `S6C_Core_Species` | 99 | Rank-biased overlap and permutation inference; Prespecified top-k composition comparisons; Shared high-loading species composition |
| S7 | `S7A_GTDB_Mapping`, `S7B_Phylo_Statistics`, `S7C_Phylo_Repro` | 342 | GTDB phylogenetic mapping; MPD, MNTD and phylogenetic signal; Cross-cohort phylogenetic reproducibility |
| S8 | `S8A_Regulatory_Resources`, `S8B_Annotation_QC` | 272 | Brain and intestinal regulatory annotation resources; Regulatory annotation technical QC |
| S9 | `S9A_Cell_GFA_Models`, `S9B_Local_Global_Align` | 248 | Complete Swedish cell-regulatory GFA model summary; Cell-local to global factor correspondence |
| S10 | `S10_HUNT_Candidates` | 22 | Complete 22-candidate HUNT validation manifest |
| S11 | `S11_HUNT_44_Directions` | 44 | Final HUNT validation across all 44 planned directions |
| S12 | `S12A_Cell_Sensitivity`, `S12B_Annotation_Overlap`, `S12C_Matched_Background_QC` | 392 | Cell-level null and sensitivity analyses; Cell-annotation overlap audit; Matched-background balance QC |
| S13 | `S13_Disease_S_LDSC_40` | 40 | Final 40 disease–annotation S-LDSC models |
| S14 | `S14A_Disease_Regions`; Supplementary Data 1: `S14B_Full_Screen` | 18,944 | Frozen disease-region definitions; comprehensive microbial regional screen supplied separately |
| S15 | `S15A_Corrected_LD`, `S15B_Coloc_Sensitivity`, `S15C_HUNT_Lookups` | 398 | Species-specific lead variants and LD pairing; Disease–species pairwise colocalization and prior sensitivity; Available HUNT locus lookups |
| S16 | `S16A_QTL_Resources`, `S16B_SCZ_TBC1D5`, `S16C_UC_EXOC3`, `S16D_CD_FUT2` | 64 | Molecular-QTL resources and input eligibility; SCZ chromosome 3 brain-cell eQTL findings; UC chromosome 5 gut-cell eQTL findings; CD–FUT2 molecular-QTL reference |
| S17 | `S17A_Variant_Harmonization`, `S17B_Common_501`, `S17C_Pairwise_Posteriors`, `S17D_Allele_Direction` | 2,268 | SCZ chromosome 3 three-way variant harmonization; Common 501-SNP comparison universe; Pairwise colocalization posteriors on the common SNP universe; Allele-direction evidence at focal and proxy variants |
| S18 | `S18A_Regional_Summary`, `S18B_Regional_Details` | 10 | Additional regional evidence summary; Additional regional molecular and sensitivity findings |

## Supplementary data attachment

Supplementary Data 1 contains the complete 18,648-row × 36-column S14B microbial regional screen. The main workbook Guide retains an S14B pointer without duplicating the full dataset. S1–S18 numbering is unchanged.

## Manuscript cross-reference audit

The frozen V3 manuscript uses the previous workbook numbering and therefore requires citation updates before submission. In particular:

- Rank-biased overlap and top-k composition results cited as former Supplementary Table 13 must be updated to **Supplementary Table S6**.
- Phylogenetic results cited as former Supplementary Table 14 must be updated to **Supplementary Table S7**.
- Former cell-model Tables 15–16 map to **S9**; former validation Tables 17–20 map to **S10–S11**; former sensitivity Tables 21–26 map to **S12**.
- The final disease S-LDSC family is now **S13** and uses the unified 40-test results; historical focused-score tables are not carried forward.
- Historical factor-GWAS/projection branches cited as former Tables 27–31 are excluded from the submission workbook and require manuscript text/citation review if retained.
- The blanket frozen-manuscript statement referring to “Supplementary Tables 1–32” is obsolete and must be changed to “Supplementary Tables S1–S18”.

The manuscript itself was not modified.

## Controlling evidence rules

- S11 uses unified BH correction across all 44 planned HUNT directions and preserves original-18 versus exploratory-four provenance.
- S13 uses the final 40-model disease results and unified q40 values while keeping original-family q values separately labelled.
- S15 preserves species-specific disease-lead/species-lead LD pairing and the updated concordant FUT2 HUNT lookup.
- S17 reports pairwise analyses only; no formal three-trait colocalization was eligible.
- S-LDSC coefficients are conditional coefficients, not marginal heritability enrichment.
