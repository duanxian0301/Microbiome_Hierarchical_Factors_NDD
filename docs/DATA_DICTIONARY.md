# Data dictionary

This repository uses tab-separated UTF-8 text tables unless otherwise noted.

## Metadata

- `dataset_accessions.tsv`: dataset, accession, DOI, genome build, sample definition, ancestry, analysis role, verification status, and provenance source.
- `software_versions.tsv`: component, recorded version, analysis scope, verification status, and provenance.
- `analysis_parameters.tsv`: final Methods parameter crosswalk. Values marked unavailable or unresolved must not be inferred.

## Integrated summary tables

- `PHASE18_HUNT_44_TEST_FINAL.tsv`: all 44 planned factor-direction records, including non-evaluable directions.
- `PHASE18_UNIFIED_22_CELL_STATUS.tsv`: status of the 22 HUNT-evaluated annotations while retaining formal and exploratory provenance.
- `PHASE18_DISEASE_FINAL_RESULTS.tsv`: the 40 final disease–annotation models.
- `PHASE18_DISEASE_UNIFIED_BH.tsv`: integrated multiplicity bookkeeping for the 40-model family.
- `PHASE18_LOCUS_MOLECULAR_EVIDENCE.tsv`: frozen regional, LD, coloc, and molecular-QTL evidence with analysis-universe distinctions.
- `PHASE18_DISEASE_COVERAGE_AUDIT.tsv`: disease coverage and source completeness audit.

## Figure source data

Subdirectories `figure1` through `figure5` contain the plotted values and annotations. They are frozen plotting inputs, not raw GWAS or QTL data. The Figure 2 files reproduce all loading, RBO, Top-k overlap, paired-rank, and factor-correspondence elements. Figure 3–5 files are the final publication-layout source tables.
