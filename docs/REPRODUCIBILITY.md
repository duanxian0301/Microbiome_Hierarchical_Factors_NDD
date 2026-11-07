# Reproducibility guide

## 1. Obtain permitted inputs

Use `metadata/dataset_accessions.tsv` to obtain each dataset directly from its provider. Do not substitute related releases without documenting the accession, genome build, and sample definition. Raw inputs are deliberately absent from this repository.

## 2. Configure paths

Set these environment variables as applicable:

- `MHFNDD_PROJECT_ROOT`: local analysis root.
- `MHFNDD_TRAIT_MANIFEST`: local trait-manifest TSV.
- `MHFNDD_COMMON_SNP_IDS`: common SNP-universe identifier file.
- `MHFNDD_HG19_TO_HG38_CHAIN`: UCSC conversion chain.
- `MHFNDD_LDSC_PYTHON`: Python 2.7 interpreter used by LDSC v1.0.1.

The public script copies preserve the executed directory relationships but contain no private machine paths. `config/paths.example.yml` is a human-readable template; historical scripts read environment variables directly.

## 3. Recorded environment

Verified software versions are listed in `metadata/software_versions.tsv`. Some historical Python and pyliftover versions were not recorded and are explicitly marked unavailable. Do not infer precise versions from a current installation.

## 4. Analysis order

The high-level execution order is:

1. harmonize cohort GWAS to the documented GRCh38, allele, ambiguity, and finite-statistic rules;
2. build the common 241-species training universe and cohort-specific error-covariance inputs;
3. fit Swedish and HUNT GFA models independently using the frozen seed and model settings;
4. align HUNT factor signs to Swedish and compare high-loading species composition;
5. construct and validate cell-regulatory annotation models;
6. run the prespecified disease S-LDSC models;
7. run targeted regional LD/coloc and molecular-QTL sensitivity analyses;
8. integrate frozen outputs and generate figure source data.

Exact script-to-result relationships are in `docs/ANALYSIS_MAP.md`. Exact frozen settings are in `metadata/analysis_parameters.tsv`.

## 5. Lightweight validation

`tests/smoke_test.py` verifies that required public files exist, key final tables have their frozen record counts, and source TSVs parse correctly. It performs no GFA, LDSC, permutation, ranking, or coloc analysis.

## 6. Limitations

The repository alone cannot recreate results that depend on non-redistributed source data. The final package is designed to make the executed transformations and frozen derived evidence inspectable after users lawfully obtain those inputs.
