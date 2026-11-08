# Microbiome Hierarchical Factors and Neuropsychiatric Disease

This repository is the reproducibility-oriented code and derived-results package for a study of host genetic organization across 241 gut microbial species in the Swedish and HUNT cohorts. The analysis combines cohort-specific global genetic factor analysis (GFA), cross-cohort comparison of high-loading species, phylogenetic organization, cell-regulatory annotations, disease stratified LD-score regression, and targeted regional and molecular-QTL follow-up.

## Scope

The repository contains the executed analysis scripts selected by the final Methods audit, small final derived result tables, and source data for the manuscript figures. It does **not** redistribute raw microbiome or disease GWAS summary statistics, LD-reference data, individual-level data, controlled-access data, or third-party molecular-QTL archives.

The Swedish and HUNT GFA models were fitted independently. Factor signs are algebraically arbitrary; HUNT factor directions were aligned to Swedish for comparison. Cross-cohort findings therefore describe recurrence of high-loading species composition and related subspace properties, not replication of fully identified biological latent entities.

## Repository layout

- `scripts/preprocessing/`: harmonization, common-universe bookkeeping, and LD-pruning scripts.
- `scripts/global_gfa/`: cohort-specific GFA fitting scripts; obtain the GFA package from its upstream source.
- `scripts/cross_cohort_validation/`: factor alignment and HUNT validation.
- `scripts/disease_association/`: cell-annotation S-LDSC and final evidence integration.
- `scripts/regional_analysis/`: targeted disease–species LD and coloc analyses.
- `scripts/figures/`: figure-data and figure build scripts.
- `metadata/`: accession, parameter, and recorded software-version crosswalks.
- `results/summary_tables/`: final integrated small result tables.
- `results/figure_source_data/`: source data for Figures 1–5.
- `docs/`: data availability, analysis map, and reproducibility guidance.
- `tests/`: lightweight integrity checks that do not rerun genome-wide analyses.

## Data access

Public accession identifiers and DOI records are listed in `metadata/dataset_accessions.tsv`. Users must obtain all primary data from the original providers and comply with their licences and access conditions. See `docs/DATA_AVAILABILITY.md` for the redistribution boundaries.

## Reproduction

The historical scripts expect the analysis directory structure described in `docs/REPRODUCIBILITY.md`. Public copies use environment variables instead of private absolute paths. Set `MHFNDD_PROJECT_ROOT` to a local analysis root after obtaining the permitted inputs. Additional path variables are required for inputs that were historically outside the project tree.

Full genome-wide reruns require substantial storage and compute. The lightweight smoke test checks package integrity and frozen table dimensions only:

```bash
python tests/smoke_test.py
```

## Scientific interpretation boundaries

- GFA factors are statistical organizations of genetic effects; they are not causal biological entities.
- M1 projection statistics are algebraic projections and are not independent latent-factor GWAS.
- Regional overlap and pairwise coloc sensitivity do not establish microbiome-mediated disease causality.
- Molecular-QTL colocalization provides cell-regulatory context and does not by itself demonstrate microbial mediation.
- Phase14A G1 was aborted and is not a completed validation result.

## Citation and licence

Author: Jing Shen. Citation metadata are provided in `CITATION.cff`; no affiliation or ORCID is asserted.

Original project code in `scripts/`, `tests/`, and `config/` is licensed under the [MIT License](LICENSE). Study-generated derived tables, figure source data, metadata, documentation, and supplementary workbooks are licensed under [Creative Commons Attribution 4.0 International](LICENSE-DATA). See [LICENSE_SCOPE.md](LICENSE_SCOPE.md) for the controlling scope and third-party exclusions.

Third-party data, code, and software are not relicensed. Raw GWAS, LD-reference, individual-level, controlled-access, and third-party molecular-QTL data are not redistributed here and remain subject to their original providers' terms.
