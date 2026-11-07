import os
from pathlib import Path
import csv, hashlib

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT = ROOT / "METHODS_FINAL_REPRODUCIBILITY_AUDIT"

def write_text(name, text):
    (OUT / name).write_text(text.strip() + "\n", encoding="utf-8")

def write_tsv(name, rows, fields):
    with (OUT / name).open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""): h.update(b)
    return h.hexdigest()

P = []
def add(section, parameter, value, source, status="DIRECTLY_VERIFIED", placement="SUPPLEMENTARY_METHODS", note=""):
    P.append(dict(section=section, parameter=parameter, exact_value=value, source=source,
                  verification_status=status, recommended_placement=placement, notes=note))

# Study design and inputs
add("Study design", "Swedish component cohorts", "SCAPIS n=8,733; SIMPLER-Västmanland n=4,515; SIMPLER-Uppsala n=981; MOS n=1,788; total n=16,017", "MANUSCRIPT_WRITING_HANDOFF/METHODS_AND_DATASETS_INVENTORY.md", placement="MAIN_METHODS")
add("Study design", "HUNT sample size", "12,652", "MANUSCRIPT_WRITING_HANDOFF/METHODS_AND_DATASETS_INVENTORY.md", placement="MAIN_METHODS")
add("Study design", "Microbial traits", "241 species shared by exact normalized species-name matching plus cohort-specific accession crosswalk; no shared sequence-level species identifier", "FULL_SCALE_GFA_FINAL_EVALUATION/inputs_241_manifest.tsv; EXTERNAL_COHORT_SPECIES_CROSSWALK.tsv", placement="MAIN_METHODS")
add("Study design", "Swedish accession range", "GCST90670368-GCST90671939 (study-wide release); 241 selected accessions listed individually in inputs_241_manifest.tsv", "MANUSCRIPT_WRITING_HANDOFF/METHODS_AND_DATASETS_INVENTORY.md; FULL_SCALE_GFA_FINAL_EVALUATION/inputs_241_manifest.tsv")
add("Study design", "HUNT accession range", "GCST90666541-GCST90667549 (study-wide release); 241 selected accessions listed individually in inputs_241_manifest.tsv", "MANUSCRIPT_WRITING_HANDOFF/METHODS_AND_DATASETS_INVENTORY.md; FULL_SCALE_GFA_FINAL_EVALUATION/inputs_241_manifest.tsv")
add("Study design", "Ancestry", "European ancestry for Swedish, HUNT, and the selected disease inputs", "METHODS_AND_DATASETS_INVENTORY.md; disease source reports", placement="MAIN_METHODS")

# Harmonization and global GFA
add("GWAS harmonization", "Analysis coordinate build", "GRCh38", "DENSE_HARMONIZED_GWAS_EXECUTION/scripts/01_harmonize_trait.py", placement="MAIN_METHODS")
add("GWAS harmonization", "Swedish coordinate conversion", "GRCh37/hg19 to GRCh38 with pyliftover and hg19ToHg38.over.chain.gz; accepted only when lifted chromosome/position exactly matched the GRCh38 reference record", "DENSE_HARMONIZED_GWAS_EXECUTION/scripts/01_harmonize_trait.py")
add("GWAS harmonization", "HUNT coordinate mode", "Direct GRCh38 in executed dense harmonization", "DENSE_HARMONIZED_GWAS_EXECUTION/scripts/01_harmonize_trait.py")
add("GWAS harmonization", "Allele matching", "Exact unordered biallelic allele-set match to the 1000G EUR GRCh38 reference; beta sign flipped when source effect allele differed from reference effect allele", "DENSE_HARMONIZED_GWAS_EXECUTION/scripts/01_harmonize_trait.py", placement="MAIN_METHODS")
add("GWAS harmonization", "Ambiguous variants", "A/T and C/G variants excluded", "DENSE_HARMONIZED_GWAS_EXECUTION/scripts/01_harmonize_trait.py; PHASE9 FULL_METHODS_AND_QC.md")
add("GWAS harmonization", "Statistic validity", "Finite beta and SE>0 required; Z=beta/SE", "DENSE_HARMONIZED_GWAS_EXECUTION/scripts/01_harmonize_trait.py")
add("GWAS harmonization", "Reference MAF filter", "1000G EUR MAF >=0.01 after LD pruning", "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA/FULL_METHODS_AND_QC.md")
add("GWAS harmonization", "Source INFO filter", "Unavailable/not applied in the project harmonization because the released species inputs did not provide a uniform INFO field", "source headers and METHODS_AND_DATASETS_INVENTORY.md", status="UNAVAILABLE")
add("GWAS harmonization", "Species sample N", "Constant cohort-level N: Swedish 16,017; HUNT 12,652; no per-SNP N", "scripts/15_fit_swedish241_rldsc.R; scripts/21_fit_hunt241_rldsc.R")
add("GWAS harmonization", "Dense common panel", "7,073,365 variants", "FULL_SCALE_GFA_FINAL_EVALUATION/common_training_identity_report.txt; Phase13 entry handoff", placement="MAIN_METHODS")
add("GWAS harmonization", "LD reference", "1000 Genomes Project Phase 3 European panel, 489 individuals, GRCh38", "PHASE9 FULL_METHODS_AND_QC.md", placement="MAIN_METHODS")
add("GWAS harmonization", "LD pruning", "PLINK2 --indep-pairwise 250kb 1 0.2", "DENSE_HARMONIZED_GWAS_EXECUTION/scripts/05_ld_prune.py", placement="MAIN_METHODS")
add("GWAS harmonization", "Training SNP count", "330,704 exact cross-cohort LD-controlled variants in Swedish dense order", "FULL_SCALE_GFA_FINAL_EVALUATION/scripts/08_build_common_training_indices.py; common_training_identity_report.txt", placement="MAIN_METHODS")
add("Error covariance", "R_ldsc input", "Cohort-specific 241x241 raw positive-definite LDSC intercept matrices; R_ldsc called with Z_hat, LD scores, M, cohort N, ncores=1, comparisons=pairs", "scripts/10_build_r241_ldsc_inputs.py; scripts/11_run_official_r241_ldsc.R")
add("Global GFA", "Executed call", "GFA::gfa_fit(Z_hat=Z, N=rep(cohort_N,241), R=R); no params argument and no F_init/fix_F/freeze_F", "scripts/15_fit_swedish241_rldsc.R; scripts/21_fit_hunt241_rldsc.R", placement="MAIN_METHODS")
add("Global GFA", "Model family", "fixed_factors error-covariance model", "patched_GFA/R/gfa_fit.R")
add("Global GFA", "Seed", "20261014", "scripts/15_fit_swedish241_rldsc.R; scripts/21_fit_hunt241_rldsc.R")
add("Global GFA", "kmax", "NULL in script/package call, resolved internally to 2*p=482", "patched_GFA/R/parameters.R; patched_GFA/R/gfa_fit.R", status="DIRECTLY_VERIFIED_PACKAGE_DEFAULT")
add("Global GFA", "Maximum iterations", "1000", "patched_GFA/R/parameters.R", status="DIRECTLY_VERIFIED_PACKAGE_DEFAULT")
add("Global GFA", "Other package defaults actually used", "cond_num=1000; max_lr_percent=0.99; lr_zero_thresh=1e-3; extrapolate=TRUE; duplicate_check_thresh=0.7; var_type=2", "patched_GFA/R/parameters.R", status="DIRECTLY_VERIFIED_PACKAGE_DEFAULT")
add("Global GFA", "Prior family", "point-normal for both F and L via flash_ebnm, optmethod='nlm'", "patched_GFA/R/parameters.R", status="DIRECTLY_VERIFIED_PACKAGE_DEFAULT")
add("Global GFA", "Initialization", "flash_greedy_init_default; fixed R factors initialized with standard-normal random matrix under the recorded seed", "patched_GFA/R/parameters.R; patched_GFA/R/gfa_fit.R", status="DIRECTLY_VERIFIED_PACKAGE_DEFAULT")
add("Global GFA", "Convergence tolerance", "Not explicitly supplied and not recorded as a scalar in the executed scripts or exported summary", "executed fit scripts and sessionInfo", status="UNAVAILABLE")
add("Global GFA", "Factor number selection", "Automatic flash greedy factor addition up to kmax, followed by null-factor removal and duplicate check; no externally chosen factor count", "patched_GFA/R/gfa_fit.R", placement="MAIN_METHODS")
add("Global GFA", "Retained global factors", "2 in Swedish and 2 in HUNT", "scripts/25_extract_and_compare_primary_factors.R; model objects", placement="MAIN_METHODS")
add("Global GFA", "Cross-cohort alignment", "Enumerated both 2-factor one-to-one assignments, selected maximum summed absolute Pearson species-loading correlation; HUNT factor signs multiplied by sign of matched Pearson correlation", "scripts/25_extract_and_compare_primary_factors.R")
add("Global GFA", "High-loading rankings", "Descending absolute loading; Top20 and Top40 taken within each factor with frozen species-order/name tie rules", "PHASE9_HUNT_VALIDATION_SELECTION_RULES.md; factor composition scripts")
add("Global GFA", "PVE", "For each trait and factor, sj/genet_var, where sj is posterior expected factor-specific standardized-effect variance and genet_var=rowSums(sj)+gvar_theta", "patched_GFA/R/pve2.R")

# Cell annotations and cell GFA
add("Cell annotations", "Brain source", "CATlas adult human brain snATAC; GEO GSE244618; 1,134,360 nuclei; 107 cell types, 42 brain regions, 3 neurotypical adult male donors; GRCh38", "PHASE9/01_CELL_ATLAS_SOURCE_MANIFEST.tsv; METHODS_AND_DATASETS_INVENTORY.md", placement="MAIN_METHODS")
add("Cell annotations", "Brain annotation definition", "Exact Table S7 cell-type cCRE assignments with Table S3 metadata; accessible chromatin cCREs, not caQTLs", "PHASE9/FULL_METHODS_AND_QC.md", placement="MAIN_METHODS")
add("Cell annotations", "Gut source", "Zhao et al. healthy-intestinal cell-type-specific DARs, Zenodo 21508422, GRCh38", "PHASE9/01_CELL_ATLAS_SOURCE_MANIFEST.tsv; METHODS_AND_DATASETS_INVENTORY.md", placement="MAIN_METHODS")
add("Cell annotations", "Gut assignment rule", "Exact published ampersand-delimited cell-token membership from Supplementary Table 4; shared memberships retained", "PHASE9/scripts/01_build_cell_atlas.py; 01_CELL_ATLAS_SOURCE_MANIFEST.tsv")
add("Cell annotations", "BED construction", "Autosomes only; converted to 0-based half-open BED; overlapping/adjacent intervals merged per cell annotation", "PHASE9/scripts/01_build_cell_atlas.py; FULL_METHODS_AND_QC.md")
add("Cell annotations", "Initial annotations", "136 total candidate cell annotations", "PHASE9/FULL_METHODS_AND_QC.md")
add("Cell annotations", "Uniform technical gate", ">=2,410 final SNPs; >=482 occupied 1-Mb bins; >=20 autosomes; no chromosome >20%; complete finite beta/SE/Z; exact ordered effect orientation; Z-matrix row rank 241", "PHASE9/FULL_METHODS_AND_QC.md")
add("Cell annotations", "Passing annotations", "124 total: 100 brain and 24 gut; 12 failed and were not rescued", "PHASE9/02_ALL_CELL_ANNOTATION_QC.tsv; FULL_METHODS_AND_QC.md", placement="MAIN_METHODS")
add("Cell annotations", "Passing SNP range", "2,475-62,084 final SNPs; median 27,852.5", "PHASE9/FULL_METHODS_AND_QC.md")
add("Cell GFA", "Fit specification", "Per annotation and cohort, same 241 species, cohort N, native positive-definite R_ldsc, seed 20261014, GFA::gfa_fit(Z_hat=Z,N=N,R=R) with unchanged package defaults", "PHASE9/scripts/05_fit_official_cell_gfa.R; FULL_METHODS_AND_QC.md", placement="MAIN_METHODS")
add("Cell GFA", "Factor selection", "Automatic package greedy/null-check procedure; retained-factor count is a model property and not itself biological support", "PHASE9/04_ALL_CELL_GFA_MODEL_SUMMARY.tsv; gfa source")
add("Cell GFA", "Global correspondence", "One-to-one assignment maximizing |cos(cell factor, global F1)|+|cos(cell factor, global F2)|; tie rules frozen; factors oriented to positive matched cosine", "PHASE9_HUNT_VALIDATION_SELECTION_RULES.md")
add("Cell GFA", "Priority score", "Within tissue: 0.35*u_R + 0.25*u_D + 0.30*u_C + 0.10*u_Q, using descending midrank percentiles", "PHASE9_HUNT_VALIDATION_SELECTION_RULES.md")

# HUNT validation
add("HUNT validation", "Original formal panel", "18 cells: six priority plus three non-priority representatives per tissue; frozen before Stage B HUNT results", "PHASE9_HUNT_VALIDATION_SELECTION_RULES.md; STAGEB_HUNT_VALIDATION_MANIFEST.tsv", placement="MAIN_METHODS")
add("HUNT validation", "Additional panel", "4 post-Stage-A Swedish-derived exploratory candidates; separate chronology and family from original 18", "STAGEB_POST_STAGEA_EXPLORATORY_MANIFEST.tsv; PHASE18_METHODS_EVIDENCE_MAP.md", placement="MAIN_METHODS")
add("HUNT validation", "Representative selection", "After excluding six priority cells per tissue, enumerate one cell per broad class; prefer all three SNP-size tertiles, otherwise maximize distinct bins; deterministic lexicographically smallest SHA256 of 20261024|region|sorted IDs", "PHASE9_HUNT_VALIDATION_SELECTION_RULES.md")
add("HUNT validation", "RBO", "Finite-list extrapolated RBO, p=0.95: (1-p)*sum[d=1..241](A_d*p^(d-1)) + A_241*p^241", "PHASE9_HUNT_VALIDATION_SELECTION_RULES.md")
add("HUNT validation", "Permutation null", "20,000 deterministic row-block permutations; seed 20261024; permutation within 9 frozen Phase7 information strata of sizes 26,23,31,29,26,25,25,31,25", "PHASE9/scripts/21_stageb_crosscohort_validation.R")
add("HUNT validation", "Permutation strata", "Phase7 primary strata derived from PC1 x PC2 tertiles of eight species information features (heritability, SE, h2 Z and R-burden measures)", "Phase7 reproducibility state; Phase9 StageB script", status="DIRECTLY_VERIFIED_DERIVED_OBJECT")
add("HUNT validation", "Empirical P", "(1 + number of null RBO >= observed RBO - 1e-15)/(20,000+1)", "PHASE9/scripts/21_stageb_crosscohort_validation.R")
add("HUNT validation", "Non-evaluable directions", "empirical P=1 in unified planned-family accounting", "Phase9 StageB script; PHASE18 build script")
add("HUNT validation", "Final multiplicity", "BH across all 44 planned factor directions (22 cells x 2); q44<=0.10", "PHASE18_HUNT_44_TEST_FINAL.tsv; PHASE18_METHODS_EVIDENCE_MAP.md", placement="MAIN_METHODS")
add("HUNT validation", "Direction support", "q44<=0.10 plus signed cosine>0 and common-Top20 sign concordance>=0.50", "PHASE18_METHODS_EVIDENCE_MAP.md", placement="MAIN_METHODS")
add("HUNT validation", "Complete support", "Both directions evaluable in both cohorts, F1/F2 preference concordant, and >=1 supported direction", "PHASE18_METHODS_EVIDENCE_MAP.md")
add("HUNT validation", "Final complete-support set", "brain__BNGA, brain__ET, brain__ITL6_2_2, brain__SMC, brain__SNCG_1, gut__Early_Colonocyte, gut__Macrophage, gut__NaiveB", "PHASE18_UNIFIED_22_CELL_STATUS.tsv", placement="MAIN_METHODS")
add("HUNT validation", "Other final statuses", "8 partial single-direction evidence; 4 not supported; 2 not evaluable", "PHASE18_UNIFIED_22_CELL_STATUS.tsv")

# Disease S-LDSC
add("Disease S-LDSC", "Diseases", "SCZ, BD, AD, CD, UC", "PHASE12A scripts; PHASE17 scripts", placement="MAIN_METHODS")
add("Disease S-LDSC", "Custom annotations", "Eight unchanged cell-associated accessible chromatin BEDs; each disease model contains exactly one custom annotation", "run_primary_phase12a.py; run_phase17.py", placement="MAIN_METHODS")
add("Disease S-LDSC", "Model", "baselineLD v2.2 plus one custom annotation; no extra generic-cell annotation and no Top1/Top5 threshold", "run_primary_phase12a.py line 164; run_phase17.py line 85", placement="MAIN_METHODS")
add("Disease S-LDSC", "Historical Phase12B Top5", "Distinct historical GLOBAL-factor Top5 analysis; not the executed configuration for the final 40 single-cell-annotation models", "PHASE12B reports; PHASE18_METHODS_EVIDENCE_MAP.md", placement="SUPPLEMENTARY_METHODS")
add("Disease S-LDSC", "LD scores", "LDSC v1.0.1; GRCh38 1000G EUR; --l2 --ld-wind-cm 1 --thin-annot; HM3/no-MHC print SNPs", "run_primary_phase12a.py; run_phase17.py")
add("Disease S-LDSC", "Regression", "--h2; --ref-ld-chr baselineLD,custom; --w-ld-chr weights.hm3_noMHC; --frqfile-chr 1000G.EUR.hg38.; --overlap-annot; --print-coefficients", "run_primary_phase12a.py; run_phase17.py", placement="MAIN_METHODS")
add("Disease S-LDSC", "Munging", "munge_sumstats.py with SNP/A1/A2, --signed-sumstats BETA,0, P, N column, and HM3 --merge-alleles", "run_primary_phase12a.py")
add("Disease S-LDSC", "Effective N", "Binary-trait Neff=4/(1/Ncase+1/Ncontrol); fixed disease counts for CD/UC; per-SNP Ncase/Ncontrol for SCZ/BD/AD where supplied", "run_primary_phase12a.py")
add("Disease S-LDSC", "Test statistic", "Custom-annotation coefficient divided by its reported coefficient SE; two-sided normal P=erfc(|z|/sqrt(2))", "parse_primary_phase12a.py; finalize_phase17.py", placement="MAIN_METHODS")
add("Disease S-LDSC", "Final family", "BH across all 40 coefficient P values as the final integrated sensitivity family; historical q30 and q10 retained", "PHASE18_DISEASE_UNIFIED_BH.tsv; PHASE18_METHODS_EVIDENCE_MAP.md", placement="MAIN_METHODS")
add("Disease S-LDSC", "Overlap-derived enrichment", "Prop._SNPs, Prop._h2 and enrichment from the historical overlap output invalidated and excluded", "PHASE12B PHASE12A_RESULT_RECONCILIATION.tsv; PHASE17 finalizer; Phase18", placement="MAIN_METHODS")
add("Disease S-LDSC", "Tau-star", "tau*sqrt[p(1-p)]*M/h2, p=annotation M_5_50/total M_5_50 and total M=5,956,110; descriptive scale only", "parse_primary_phase12a.py; finalize_phase17.py", placement="SUPPLEMENTARY_METHODS")

# Regional and molecular-QTL analyses
add("Regional analysis", "Disease-risk regions", "296 fixed regions: SCZ 130, AD 41, CD 77, UC 48", "PHASE13_ENTRY_HANDOFF.md; Phase13 execution log", placement="MAIN_METHODS")
add("Regional analysis", "Region derivation", "Disease-only HM3 clumping P1=5e-8, r2=0.1, 1-Mb window; lead +/-250 kb; merge overlap or gaps <=100 kb", "PHASE13/scripts/run_phase13b.py", placement="SUPPLEMENTARY_METHODS")
add("Regional analysis", "Species sets", "Primary Swedish F1/F2 Top20 union=31 species; sensitivity Top40 union=63 species", "PHASE13B_FOUR_DISEASE_LOCUS_REPORT.md", placement="MAIN_METHODS")
add("Regional analysis", "Screening rows", "296 x 63=18,648 region-species combinations; disease regions fixed before species P inspection", "PHASE13 execution log; Phase13B report")
add("Regional analysis", "Same-LD support", "Species regional P<1e-4 and EUR r2>=0.1 between that species' own regional lead and the disease lead", "PHASE14B_DISEASE_REGION_INTEGRATION.md", placement="MAIN_METHODS")
add("Regional analysis", "LD for local pairing", "1000G Phase 3 EUR, 489 samples; species-specific paired r2", "PHASE14/scripts/run_phase14b.py")
add("Regional analysis", "Six fixed targets", "SCZ chr3 and chr11; CD chr6 and FUT2 chr19; UC chr5 and chr1; fixed by Phase14 protocol, no P-driven expansion", "PHASE14_ENTRY_HANDOFF.md; Phase14 reports", placement="MAIN_METHODS")
add("Colocalization", "Method/version", "coloc.abf, coloc 6.0.1", "PHASE16_ANALYSIS_QC.md", placement="MAIN_METHODS")
add("Colocalization", "Default priors", "p1=1e-4, p2=1e-4, p12=1e-5", "Phase14 run_coloc_abf.R defaults; Phase16 explicit script", placement="MAIN_METHODS")
add("Colocalization", "Sensitivity", "p12=1e-6, 1e-5, 1e-4 with p1=p2=1e-4", "Phase14/16 coloc scripts")
add("Colocalization", "Minimum regional SNPs", ">=100 harmonized common SNPs for Phase14 disease-species coloc", "Phase14 run_coloc_abf.R")
add("Colocalization", "Species input limitations", "1000G EUR reference MAF, cohort-level N=16,017; original species INFO and SNP-specific N unavailable", "Phase14 run_coloc_abf.R; Phase16 ENTRY_AUDIT.md", placement="MAIN_METHODS")
add("Colocalization", "SCZ chr3 Phase14 universe", "1,030 harmonized SCZ-species SNPs after reference-MAF filtering", "Phase16 ENTRY_AUDIT.md")
add("Colocalization", "SCZ chr3 three-pair common universe", "501 shared non-palindromic SNPs per Bryois cell type after excluding 111 ambiguous variants from 612 overlaps", "Phase16 ENTRY_AUDIT.md", placement="MAIN_METHODS")
add("Molecular QTL", "Bryois", "Zenodo 7276971; 8 brain cell types; n=144-192; hg19 and GRCh38 positions; beta and nominal P; SE reconstructed from beta/P; no MAF/INFO/per-SNP N/matched LD", "Phase15 CELL_QTL_RESOURCE_AUDIT.tsv; Phase16 ENTRY_AUDIT.md", placement="MAIN_METHODS")
add("Molecular QTL", "CEDAR-2", "Zenodo 21457760; 57 healthy analyzed; 401 hierarchy nodes; GRCh38; conditional and permutational beta/P/SE/q/alleles; no matched LD", "Phase15 CELL_QTL_RESOURCE_AUDIT.tsv", placement="MAIN_METHODS")
add("Molecular QTL", "Allele harmonization", "Match chromosome/position and allele set; align beta to a common effect allele; exclude unresolved palindromic variants when source-specific allele frequencies absent", "Phase16 phase16_prepare.py; ENTRY_AUDIT.md")
add("Molecular QTL", "Formal three-trait analysis", "Not performed: HyPrColoc/moloc unavailable and study-specific MAF/INFO/per-SNP N, matched LD, and single-signal assumptions were not adequate", "Phase16 PHASE16_FINAL_SCIENTIFIC_REPORT.md", placement="MAIN_METHODS")

# Multiple testing / availability
add("Multiplicity", "Global HUNT cell validation", "BH q<=0.10 across 44 planned directions", "Phase18")
add("Multiplicity", "Disease S-LDSC", "BH q<0.05 across final 40 coefficient tests; q30/q10 preserved historically", "Phase18")
add("Multiplicity", "Phase13 regional background", "Empirical background was not a formal no-sharing null; all matched-background q=1 and not used as confirmatory evidence", "PHASE14_ENTRY_HANDOFF.md")
add("Reproducibility", "Historical-output policy", "All analyses read from frozen scripts, logs, manifests and objects; this audit reran no statistical analysis", "current audit protocol", placement="MAIN_METHODS")

write_tsv("METHODS_PARAMETER_CROSSWALK.tsv", P,
          ["section","parameter","exact_value","source","verification_status","recommended_placement","notes"])

datasets = [
dict(dataset="Swedish microbiome GWAS", accession="GCST90670368-GCST90671939", doi="10.1038/s41588-026-02512-2", build="GRCh37 source; lifted to GRCh38", sample="16,017", ancestry="European", role="Discovery microbiome cohort", status="DIRECTLY_VERIFIED", source="METHODS_AND_DATASETS_INVENTORY.md; inputs_241_manifest.tsv"),
dict(dataset="HUNT microbiome GWAS", accession="GCST90666541-GCST90667549", doi="10.1038/s41588-026-02502-4", build="GRCh38 in project input (release provenance also reports GRCh37-to-GRCh38 standardization)", sample="12,652", ancestry="European", role="Cross-cohort composition validation", status="DISCREPANCY_FLAGGED", source="01_harmonize_trait.py; METHODS_AND_DATASETS_INVENTORY.md"),
dict(dataset="SCZ PGC3 wave3", accession="PGC3_SCZ_wave3.european.autosome.public.v3; Figshare 19426775", doi="10.1038/s41586-022-04434-5", build="GRCh37 source; harmonized GRCh38", sample="52,017 cases; 75,889 controls; 1,369 trios (per-variant maxima 53,386/77,258)", ancestry="European", role="Disease GWAS", status="DIRECTLY_VERIFIED", source="DISEASE_SOURCE_ACCESS_AUDIT.tsv; METHODS_AND_DATASETS_INVENTORY.md"),
dict(dataset="BD PGC4", accession="bip2024_eur_no23andMe; Figshare 27216117", doi="10.1038/s41586-024-08468-9", build="GRCh37 source; harmonized GRCh38", sample="59,287 cases; 781,022 controls", ancestry="European", role="Disease GWAS", status="DIRECTLY_VERIFIED", source="DISEASE_SOURCE_ACCESS_AUDIT.tsv; METHODS_AND_DATASETS_INVENTORY.md"),
dict(dataset="AD EADB/IGAP/PGC no-biobank", accession="GCST90704648", doi="10.1038/s41588-026-02583-1", build="GRCh38", sample="Case/control totals not verified in local audit", ancestry="European", role="Disease GWAS", status="PARTIALLY_VERIFIED", source="DISEASE_SOURCE_ACCESS_AUDIT.tsv; METHODS_AND_DATASETS_INVENTORY.md"),
dict(dataset="Crohn disease", accession="GCST004132", doi="10.1038/ng.3760", build="GRCh37 source; harmonized GRCh38", sample="12,194 cases; 28,072 controls", ancestry="European", role="Disease GWAS", status="DIRECTLY_VERIFIED", source="DISEASE_SOURCE_ACCESS_AUDIT.tsv"),
dict(dataset="Ulcerative colitis", accession="GCST004133", doi="10.1038/ng.3760", build="GRCh37 source; harmonized GRCh38", sample="12,366 cases; 33,609 controls", ancestry="European", role="Disease GWAS", status="DIRECTLY_VERIFIED", source="DISEASE_SOURCE_ACCESS_AUDIT.tsv"),
dict(dataset="1000 Genomes Phase 3 EUR", accession="1000G EUR GRCh38 PLINK bundle", doi="", build="GRCh38", sample="489", ancestry="European", role="Allele/MAF/LD reference", status="DIRECTLY_VERIFIED", source="PHASE9 FULL_METHODS_AND_QC.md"),
dict(dataset="baselineLD", accession="baselineLD v2.2", doi="", build="GRCh38", sample="Reference annotation bundle", ancestry="", role="S-LDSC baseline annotations", status="DIRECTLY_VERIFIED", source="Phase12A/17 scripts"),
dict(dataset="CATlas adult brain snATAC", accession="GEO GSE244618", doi="10.1126/science.adf7044", build="GRCh38", sample="1,134,360 nuclei; 3 neurotypical adult male donors", ancestry="Not used for GWAS", role="Cell-associated accessible cCRE annotations", status="DIRECTLY_VERIFIED", source="PHASE9 source manifest; METHODS_AND_DATASETS_INVENTORY.md"),
dict(dataset="Zhao healthy intestine atlas", accession="Zenodo 21508422", doi="10.1038/s41588-026-02755-z", build="GRCh38", sample="Cell/donor counts annotation-specific; healthy controls", ancestry="Not used for GWAS", role="Cell-type DAR annotations", status="DIRECTLY_VERIFIED", source="PHASE9 source manifest"),
dict(dataset="Bryois brain cell eQTL", accession="Zenodo 7276971; EGA EGAS00001006345", doi="10.1038/s41593-022-01128-z", build="hg19 and GRCh38 supplied", sample="144-192 by cell type", ancestry="", role="SCZ chr3 molecular-QTL context", status="DIRECTLY_VERIFIED", source="Phase15 CELL_QTL_RESOURCE_AUDIT.tsv"),
dict(dataset="CEDAR-2 gut cell eQTL", accession="Zenodo 21457760", doi="10.1038/s41467-026-76672-4", build="GRCh38", sample="57 healthy analyzed (60 collected)", ancestry="", role="UC chr5 and FUT2 molecular-QTL context", status="DIRECTLY_VERIFIED", source="Phase15 CELL_QTL_RESOURCE_AUDIT.tsv"),
dict(dataset="GENCODE", accession="GENCODE v48 basic GTF", doi="", build="GRCh38", sample="NA", ancestry="", role="Gene annotation", status="DIRECTLY_VERIFIED", source="Phase15 CELL_QTL_RESOURCE_AUDIT.tsv"),
]
write_tsv("METHODS_DATASET_ACCESSION_CROSSWALK.tsv", datasets,
          ["dataset","accession","doi","build","sample","ancestry","role","status","source"])

software = [
dict(component="R", version="4.4.0 (2024-04-24 ucrt)", scope="Global/cell GFA and R_ldsc", status="DIRECTLY_VERIFIED", source="Swedish model sessionInfo; R_ldsc sessionInfo"),
dict(component="GFA", version="1.0.0.0470 as sessionInfo; narrative tables abbreviate 1.0.0.470", scope="Global and cell GFA", status="DIRECTLY_VERIFIED_WITH_DISPLAY_DISCREPANCY", source="model sessionInfo; PHASE9 model summary"),
dict(component="GFA source commit", version="eb582bb28725b28dd8b4b188da8579288a9e69b7", scope="Cell GFA", status="DIRECTLY_VERIFIED", source="PHASE9 _stagea_model_audit_checkpoint.tsv"),
dict(component="flashier", version="1.0.60", scope="GFA", status="DIRECTLY_VERIFIED", source="model sessionInfo"),
dict(component="ebnm", version="1.1-42", scope="GFA", status="DIRECTLY_VERIFIED", source="model sessionInfo"),
dict(component="PLINK2", version="v2.0.0-b.1-dev 64-bit (7 Oct 2026)", scope="LD pruning", status="DIRECTLY_VERIFIED", source="FULL_SCALE_GFA_FINAL_EVALUATION/ld_prune/*.log"),
dict(component="LDSC", version="v1.0.1", scope="Disease S-LDSC", status="DIRECTLY_VERIFIED", source="Phase12A/17 logs and methods"),
dict(component="Python runtime for LDSC", version="Python 2.7 bundled interpreter", scope="LDSC/munge_sumstats", status="DIRECTLY_VERIFIED_PATH_ONLY", source="run_primary_phase12a.py PY27 path"),
dict(component="Python runtime for project orchestration", version="Not recorded consistently in historical logs", scope="Harmonization/audit scripts", status="UNAVAILABLE", source="executed Python scripts/logs"),
dict(component="pyliftover", version="Not recorded", scope="GRCh37-to-GRCh38 conversion", status="UNAVAILABLE", source="01_harmonize_trait.py"),
dict(component="UCSC hg19ToHg38 chain", version="local file hg19ToHg38.over.chain.gz; checksum should be cited from input manifest", scope="Coordinate conversion", status="DIRECTLY_VERIFIED_FILE", source="01_harmonize_trait.py"),
dict(component="coloc", version="6.0.1", scope="Regional and molecular-QTL coloc.abf", status="DIRECTLY_VERIFIED", source="PHASE16_ANALYSIS_QC.md"),
dict(component="susieR", version="0.14.2", scope="Eligibility audit only; SuSiE not run", status="DIRECTLY_VERIFIED", source="PHASE16_ANALYSIS_QC.md"),
dict(component="HyPrColoc", version="Not installed", scope="Three-trait eligibility audit", status="DIRECTLY_VERIFIED", source="PHASE16_ANALYSIS_QC.md"),
dict(component="moloc", version="Not installed", scope="Three-trait eligibility audit", status="DIRECTLY_VERIFIED", source="PHASE16_ANALYSIS_QC.md"),
]
write_tsv("METHODS_SOFTWARE_VERSION_CROSSWALK.tsv", software,
          ["component","version","scope","status","source"])

main = r'''# Main Methods specification (structured; not manuscript prose)

## Scope

This file is a claim-to-parameter specification for drafting. It is not manuscript prose and does not replace the row-level crosswalk.

## 1. Cohorts and microbial traits

- State the four Swedish component cohorts and n=16,017; state HUNT n=12,652.
- Define the 241-species crosswalk as exact normalized names plus cohort-specific accessions; do not imply sequence-ID identity.
- Identify Swedish as discovery and HUNT as cross-cohort composition validation.

## 2. Harmonization and global GFA

- State GRCh38 analysis coordinates, exact allele-set matching, effect-sign flipping, ambiguous-variant exclusion, 1000G EUR MAF>=0.01, and the absence of a uniform source INFO filter.
- State 7,073,365 dense and 330,704 LD-controlled training variants, 1000G EUR n=489, and PLINK2 250 kb/step 1/r2=0.2 pruning.
- Specify the executed GFA call, cohort-specific R_ldsc, seed 20261014, automatic factor selection, and two retained factors in each cohort.
- Do not describe package defaults as user-tuned parameters. The exact defaults belong in Supplementary Methods.

## 3. Cell-regulatory GFA and HUNT validation

- Distinguish CATlas accessible cCREs and Zhao DARs from caQTLs.
- State that 124/136 annotations passed the uniform technical gate (100 brain, 24 gut) and were fitted with the same GFA framework.
- Preserve chronology: original 18-cell panel first; four post-Stage-A exploratory cells second.
- State RBO p=0.95, 20,000 stratified permutations, seed 20261024, unified BH across 44 planned directions, and the complete-support rule.
- Name the final eight complete-support annotations. Do not infer biology from factor retention alone.

## 4. Disease S-LDSC

- Specify exactly baselineLD v2.2 plus one unchanged custom cell annotation per model, HM3/no-MHC, 1000G EUR GRCh38, 1-cM LD-score windows, and coefficient-based inference.
- State that the final integrated BH family contains 40 coefficient tests. Preserve historical q30/q10 only as provenance.
- Exclude overlap-derived Prop_h2/enrichment. Tau-star is descriptive and belongs in Supplementary Methods.
- Do not describe the distinct Phase12B GLOBAL Top5 analysis as the final 40-model design.

## 5. Regional and molecular-QTL analysis

- State the fixed 296-region derivation and the Top20-union primary/Top40-union sensitivity species sets.
- Require species-specific lead-to-disease-lead LD pairing; do not transfer an r2 or coloc result between species.
- Describe the six Phase14 targets as fixed follow-up regions, including FUT2 as a known reference.
- Give coloc.abf default priors and p12 sensitivity; disclose reference MAF, study-level species N, missing INFO/SNP-N, and the single-signal limitation.
- Distinguish disease-species, disease-eQTL, and species-eQTL comparisons. State explicitly that no formal three-trait posterior was calculated.

## 6. Reproducibility statement

- Cite the frozen scripts, logs, model objects, SHA-256 manifests, and version crosswalk.
- Describe post hoc/exploratory families explicitly: the four extra HUNT cells, integrated q44, integrated S-LDSC q40, and local/molecular-QTL follow-up.
'''
write_text("METHODS_MAIN_TEXT_SPECIFICATION.md", main)

supp = r'''# Supplementary Methods specification (structured; not manuscript prose)

## A. Species GWAS harmonization

Report source column mappings, liftOver chain path/hash, coordinate acceptance rule, unordered allele-set comparison, beta flipping, non-finite-statistic exclusion, ambiguous-variant exclusion, reference MAF threshold, and the lack of uniform INFO/per-SNP N. Include dense/training manifest hashes and per-cohort accession lists.

## B. R_ldsc and GFA

Record the exact R_ldsc call and matrix QC. Reproduce all actual GFA defaults from `parameters.R`: kmax NULL->482, cond_num 1000, max_lr_percent 0.99, lr_zero_thresh 1e-3, max_iter 1000, extrapolate TRUE, point-normal EBNM with nlm for F/L, greedy default initialization, duplicate threshold 0.7, var_type 2. State that convergence tolerance was not explicitly recorded. Give the PVE formula and cross-cohort factor permutation/sign rule.

## C. Cell annotation construction

For CATlas, document Table S7 cCRE assignments and Table S3 metadata. For Zhao, document healthy-only Supplementary Table 4 DARs and exact multi-label token splitting. Give coordinate conventions, autosome filtering, interval merging, all uniform gate thresholds, per-cell SNP/QC/M_5_50 values, and BED hashes.

## D. Cell factor correspondence and HUNT testing

Give the cosine assignment/tie/orientation rules, priority-score components and weights, representative-selection hash rule, exact RBO formula, Phase7 strata sizes, empirical-P formula, non-evaluable P=1 treatment, unified q44 calculation, signed cosine/common-Top20 criterion, preference-concordance rule, and all 22 final cell statuses. Keep the original 18 and exploratory four in separate columns.

## E. Disease S-LDSC command specification

Include the exact `--l2` and `--h2` argument sets, annotation builder, HM3/no-MHC lists, baseline/weights/frequency paths, disease-specific N construction, coefficient z/P calculation, q30/q10 provenance, q40 integrated sensitivity, corrected M_5_50 tau-star formula, annotation QC and collinearity note. Document why the overlap-derived enrichment fields were invalidated.

## F. Regional analyses

Document disease-only clumping (P1=5e-8, r2=0.1, 1 Mb), lead +/-250 kb regions, merge gap <=100 kb, disease counts, 31/63 species unions, regional P<1e-4 plus paired r2>=0.1 rule, 1000G EUR LD, all six fixed targets, and the Phase13 pairing defect corrected in Phase14.

## G. Colocalization and molecular-QTL

Report coloc 6.0.1, p1/p2/p12, minimum 100 SNP Phase14 eligibility, exact per-comparison SNP universes, allele harmonization, reference-MAF substitution, sample-N limitations, Bryois SE reconstruction, CEDAR-2 conditional/permutational input distinctions, and absence of matched LD. State that pairwise PP4 cannot be combined into a three-trait posterior and SuSiE/HyPrColoc/moloc were not qualified.

## H. Availability and audit trail

List public study accessions/DOIs, controlled-access resources, locally derived but not yet publicly archived code/data, input/output SHA-256 files, and every unresolved item. No claim should imply that local derived products already have a permanent public repository accession.
'''
write_text("METHODS_SUPPLEMENTARY_SPECIFICATION.md", supp)

unresolved = r'''# Unresolved and discrepancy items

## Must be resolved or disclosed before submission

1. **AD sample totals:** the selected GCST90704648 file and build are verified, but exact final case/control totals were not verified locally. Do not insert totals from a related meta-analysis without source-file confirmation.
2. **HUNT source-build wording:** the project inventory describes source GRCh37 followed by GRCh38 standardization, whereas the executed dense harmonizer treats the selected HUNT files as direct GRCh38. The final Methods should describe the actual input file/build after checking the release metadata, and retain the discrepancy note until then.
3. **GFA version formatting:** actual `sessionInfo()` reports `GFA_1.0.0.0470`; several reports abbreviate it to `1.0.0.470`. Use the sessionInfo form or state both as the same local package build.
4. **GFA convergence tolerance:** `max_iter=1000` and `extrapolate=TRUE` are verified, but no scalar stopping tolerance was explicitly passed or exported. Do not infer one from a current flashier installation.
5. **Python/pyliftover versions:** historical orchestration logs do not consistently record Python or pyliftover versions. Paths and script hashes are available; version numbers are unavailable.
6. **Chain-file checksum:** the exact chain path is verified. Ensure its SHA-256 is added to the public reproducibility bundle if the file can be redistributed.
7. **Species GWAS metadata:** per-SNP N, INFO and study-specific MAF are not uniformly available. Regional coloc used cohort-level N and 1000G EUR reference MAF and must remain sensitivity-grade.
8. **Disease INFO:** CD/UC INFO is unavailable in the regional inputs; other disease files have heterogeneous metadata.
9. **Phase12A overlap-derived fields:** historical Prop._SNPs/Prop._h2/enrichment output was invalidated by the re-audit and must not reappear. The coefficient, SE and coefficient P are the inferential objects.
10. **Phase12B description conflict:** Phase12B GLOBAL Top5 analyses are distinct from the final 40 single-cell-annotation models. Any historical text that says the final cell disease models used Top5 or a generic-cell covariate must be corrected in future manuscript drafting, without altering historical files.
11. **Multiplicity chronology:** q44 and q40 are later integrated harmonizations. They are the recommended final families, but their post hoc consolidation must be disclosed; original q36/q8 and q30/q10 remain provenance.
12. **Molecular-QTL LD and QC:** no matched study LD, complete per-variant INFO, or per-variant N is available for the focal species/Bryois/CEDAR-2 comparisons. Multi-signal SuSiE and formal three-trait colocalization were therefore not qualified.
13. **CEDAR-2 archive scope:** only the necessary conditional and permutational archives were downloaded; the full nominal archive was not. State this when defining the local SNP universe.
14. **Public archiving:** many derived tables, scripts and model objects are local and SHA-256-tracked but do not yet have permanent public repository accessions. A code/data release plan is still required.
15. **Factor polarity:** factor signs are algebraically arbitrary. HUNT was aligned to Swedish by loading correlation; no biological positive-direction anchor was imposed on Swedish. Use absolute loadings for membership and explicit signed loadings for direction.
16. **Phase14A G1:** the leakage-free refit was aborted and is not a scientific result. Candidate 54 and incomplete artifacts must not enter Methods/Results as completed validation.
'''
write_text("METHODS_UNRESOLVED_ITEMS.md", unresolved)

input_paths = [
ROOT/"MANUSCRIPT_WRITING_HANDOFF"/"METHODS_AND_DATASETS_INVENTORY.md",
ROOT/"MANUSCRIPT_INTEGRATION_FREEZE_20261031"/"REVISION_3"/"MASTER_RESEARCH_HANDOFF_V3.md",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"scripts"/"08_build_common_training_indices.py",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"scripts"/"10_build_r241_ldsc_inputs.py",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"scripts"/"11_run_official_r241_ldsc.R",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"scripts"/"15_fit_swedish241_rldsc.R",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"scripts"/"21_fit_hunt241_rldsc.R",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"scripts"/"25_extract_and_compare_primary_factors.R",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"patched_GFA"/"R"/"parameters.R",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"patched_GFA"/"R"/"gfa_fit.R",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"patched_GFA"/"R"/"pve2.R",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"models"/"Swedish241_genomewide_ldcontrolled_Rldsc_seed20261014_sessionInfo.txt",
ROOT/"CELL_ANNOTATION_FACTOR_LOCALIZATION"/"DENSE_HARMONIZED_GWAS_EXECUTION"/"scripts"/"01_harmonize_trait.py",
ROOT/"CELL_ANNOTATION_FACTOR_LOCALIZATION"/"DENSE_HARMONIZED_GWAS_EXECUTION"/"scripts"/"05_ld_prune.py",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"01_CELL_ATLAS_SOURCE_MANIFEST.tsv",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"02_ALL_CELL_ANNOTATION_QC.tsv",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"04_ALL_CELL_GFA_MODEL_SUMMARY.tsv",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"PHASE9_HUNT_VALIDATION_SELECTION_RULES.md",
ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"scripts"/"21_stageb_crosscohort_validation.R",
ROOT/"PHASE12A_CELL_DISEASE_SLD SC"/"scripts"/"run_primary_phase12a.py",
ROOT/"PHASE12A_CELL_DISEASE_SLD SC"/"scripts"/"parse_primary_phase12a.py",
ROOT/"PHASE17_EXPLORATORY_CELL_DISEASE_EXTENSION"/"scripts"/"run_phase17.py",
ROOT/"PHASE17_EXPLORATORY_CELL_DISEASE_EXTENSION"/"scripts"/"finalize_phase17.py",
ROOT/"PHASE18_MANUSCRIPT_EVIDENCE_FINAL"/"PHASE18_HUNT_44_TEST_FINAL.tsv",
ROOT/"PHASE18_MANUSCRIPT_EVIDENCE_FINAL"/"PHASE18_UNIFIED_22_CELL_STATUS.tsv",
ROOT/"PHASE18_MANUSCRIPT_EVIDENCE_FINAL"/"PHASE18_DISEASE_FINAL_RESULTS.tsv",
ROOT/"PHASE13_SHARED_STRUCTURE_AND_DISEASE_LOCI"/"scripts"/"run_phase13b.py",
ROOT/"PHASE14_TARGETED_VALIDATION"/"B_SCZ_CD_UC_LOCI"/"scripts"/"run_phase14b.py",
ROOT/"PHASE14_TARGETED_VALIDATION"/"B_SCZ_CD_UC_LOCI"/"scripts"/"run_coloc_abf.R",
ROOT/"PHASE15_CELL_MOLECULAR_AND_MANUSCRIPT"/"CELL_QTL_RESOURCE_AUDIT.tsv",
ROOT/"PHASE16_SCZ_THREE_WAY_SIGNAL_AUDIT"/"ENTRY_AUDIT.md",
ROOT/"PHASE16_SCZ_THREE_WAY_SIGNAL_AUDIT"/"scripts"/"run_phase16_pairwise_coloc.R",
ROOT/"PHASE16_SCZ_THREE_WAY_SIGNAL_AUDIT"/"PHASE16_ANALYSIS_QC.md",
]
manifest=[]
for p in input_paths:
    manifest.append(dict(path=str(p), exists=str(p.exists()).upper(), bytes=p.stat().st_size if p.exists() else "NA",
                         sha256=sha(p) if p.exists() else "NA", role="audited source"))
write_tsv("INPUT_SHA256_MANIFEST.tsv", manifest, ["path","exists","bytes","sha256","role"])

out_names=["METHODS_MAIN_TEXT_SPECIFICATION.md","METHODS_SUPPLEMENTARY_SPECIFICATION.md",
"METHODS_PARAMETER_CROSSWALK.tsv","METHODS_UNRESOLVED_ITEMS.md",
"METHODS_DATASET_ACCESSION_CROSSWALK.tsv","METHODS_SOFTWARE_VERSION_CROSSWALK.tsv",
"INPUT_SHA256_MANIFEST.tsv","scripts/build_methods_audit.py"]
rows=[]
for name in out_names:
    p=OUT/name; rows.append(dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)))
write_tsv("OUTPUT_SHA256_MANIFEST.tsv", rows, ["path","bytes","sha256"])

print(f"wrote {len(out_names)+1} deliverables; parameters={len(P)}; inputs={len(manifest)}")
