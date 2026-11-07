#!/usr/bin/env python3
"""Parse, QC, multiplicity-correct, document, and hash Phase17 outputs."""
from __future__ import annotations
import csv, gzip, hashlib, math, re
from collections import Counter
import os
from pathlib import Path

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT=ROOT/"PHASE17_EXPLORATORY_CELL_DISEASE_EXTENSION"
P9=ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"
P12=ROOT/"PHASE12A_CELL_DISEASE_SLD SC"
P12B=ROOT/"PHASE12B_MICROBIOME_INFORMED_DISEASE_INTEGRATION"
REF=ROOT/"GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT"/"CELL_ANNOTATION_C1_EXECUTION"/"reference_download"/"GRCh38"/"baselineLD_v2.2_hm3"
CELLS=["brain__ET","brain__ITL6_2_2"]
DISEASES=["SCZ","BD","AD","CD","UC"]
FORMAL=["brain__BNGA","brain__SNCG_1","brain__SMC","gut__Macrophage","gut__Early_Colonocyte","gut__NaiveB"]

def read_tsv(p):
  with open(p,encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h,delimiter="\t"))
def write_tsv(p,rows,fields=None):
  if fields is None:fields=list(rows[0])
  with open(p,"w",encoding="utf-8",newline="") as h:
    w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",extrasaction="ignore",lineterminator="\n");w.writeheader();w.writerows(rows)
def first(p):return float(p.read_text().split()[0])
def bh(vals):
  o=sorted(range(len(vals)),key=lambda i:vals[i]);q=[1.0]*len(vals);run=1.0;n=len(vals)
  for rank in range(n,0,-1):
    i=o[rank-1];run=min(run,vals[i]*n/rank);q[i]=run
  return q
def fmt(x):return "%.10g"%x
def sha(p):
  h=hashlib.sha256()
  with open(p,"rb") as f:
    for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
  return h.hexdigest()

def parse_log(path):
  s=path.read_text(errors="replace")
  h2,se=map(float,re.findall(r"Total Observed scale h2: ([0-9.eE+\-]+) \(([0-9.eE+\-]+)\)",s)[-1])
  stem=path.name.removesuffix(".h2.log")
  hit=[x for x in read_tsv(OUT/"results"/(stem+".results")) if x["Category"]=="L2_1"]
  if len(hit)!=1:raise RuntimeError("expected one L2_1 row in "+stem)
  c=float(hit[0]["Coefficient"]);e=float(hit[0]["Coefficient_std_error"]);z=c/e;p=math.erfc(abs(z)/math.sqrt(2))
  return c,e,z,p,h2,se,("YES" if "ill-conditioned" in s.lower() else "NO"),s

def inventory():
  models=read_tsv(P9/"04_ALL_CELL_GFA_MODEL_SUMMARY.tsv")
  formal={x["annotation_id"] for x in read_tsv(P9/"STAGEB_HUNT_VALIDATION_MANIFEST.tsv")}
  exp={x["annotation_id"]:x for x in read_tsv(P9/"STAGEB_POST_STAGEA_EXPLORATORY_RESULTS.tsv")}
  rows=[]
  for x in models:
    cid=x["annotation_id"]
    rows.append({"annotation_id":cid,"organ":x["organ"],"model_status":x["model_status"],"retained_factors":x["retained_factors"],"technical_valid":x["technical_valid"],"phase9_formal_hunt":str(cid in formal).upper(),"phase9_exploratory_hunt":str(cid in exp).upper(),"phase9_cross_cohort_support":("FORMAL_SUPPORTED" if cid in FORMAL else (exp[cid]["validation_status"] if cid in exp else "NO_FORMAL_OR_EXPLORATORY_SUPPORT")),"phase12a_disease_sldsc":str(cid in FORMAL).upper(),"phase17_target":str(cid in CELLS).upper()})
  write_tsv(OUT/"PHASE17_CELL_INVENTORY.tsv",rows)
  return models,formal,exp

def main():
  models,formal,exp=inventory();meta={x["annotation_id"]:x for x in read_tsv(P9/"01_CELL_ATLAS_SOURCE_MANIFEST.tsv")}
  mtotal=sum(first(REF/("baselineLD.%d.l2.M_5_50"%c)) for c in range(1,23))
  qc=[];sets={}
  for cell in CELLS:
    aset=set();total=0;allm=m550=0;complete=True
    for c in range(1,23):
      pre=OUT/"ldscores"/(cell+".%d"%c)
      complete &= all(Path(str(pre)+s).exists() for s in [".l2.ldscore.gz",".l2.M",".l2.M_5_50"])
      allm+=first(Path(str(pre)+".l2.M"));m550+=first(Path(str(pre)+".l2.M_5_50"))
      with gzip.open(OUT/"ldscores"/(cell+".%d.annot.gz"%c),"rt") as h:
        next(h)
        for z in h:
          total+=1;f=z.rstrip().split("\t")
          if f[-1]=="1":aset.add(f[2])
    sets[cell]=aset;r=meta[cell]
    qc.append({"CELL_ID":cell,"COORDINATE_BUILD":r["coordinate_build"],"SOURCE_DEFINITION":r["annotation_definition"],"SOURCE_BED":r["bed_path"],"SOURCE_BED_SHA256_EXPECTED":r["bed_sha256"],"SOURCE_BED_SHA256_OBSERVED":sha(Path(r["bed_path"])),"ORIGINAL_ASSIGNMENT_PEAKS":r["original_assignment_peaks"],"MERGED_AUTOSOMAL_INTERVALS":r["merged_autosomal_intervals"],"MERGED_AUTOSOMAL_BP":r["merged_autosomal_bp"],"HM3_BASELINE_ROWS":total,"HM3_ANNOTATED":len(aset),"HM3_COVERAGE":fmt(len(aset)/total),"M_ALL":fmt(allm),"M_5_50":fmt(m550),"M_5_50_PROP":fmt(m550/mtotal),"LDSCORE_CHR_COMPLETE":str(complete).upper(),"QC_STATUS":"PASS" if complete and sha(Path(r["bed_path"]))==r["bed_sha256"] and len(aset)>0 else "FAIL"})
  inter=len(sets[CELLS[0]]&sets[CELLS[1]]);union=len(sets[CELLS[0]]|sets[CELLS[1]])
  for r in qc:r.update(TARGET_PAIR_HM3_INTERSECTION=inter,TARGET_PAIR_JACCARD=fmt(inter/union),COLLINEARITY_NOTE="Pair overlap reported; each disease model contains one custom annotation plus baselineLD, not both target annotations jointly.")
  write_tsv(OUT/"PHASE17_ANNOTATION_QC.tsv",qc)
  if any(x["QC_STATUS"]!="PASS" for x in qc):raise RuntimeError("annotation QC failed")
  new=[]
  for d in DISEASES:
    for cell in CELLS:
      c,se,z,p,h2,h2se,ill,text=parse_log(OUT/"logs"/(d+"__"+cell+".h2.log"))
      ma=float(next(x["M_5_50"] for x in qc if x["CELL_ID"]==cell));prop=ma/mtotal;sd=math.sqrt(prop*(1-prop))
      new.append({"FAMILY":"ADDITIONAL_EXPLORATORY_10","DISEASE":d,"CELL_ID":cell,"CELL_LABEL":meta[cell]["cell_type_label"],"ANNOTATION_CLASS":"CELL_ASSOCIATED_ACCESSIBLE_CCRE","MODEL":"baselineLD_v2.2+single_custom_annotation","STATUS":"COMPLETED","H2_OBS":fmt(h2),"H2_OBS_SE":fmt(h2se),"COEFFICIENT":fmt(c),"COEFFICIENT_SE":fmt(se),"COEFFICIENT_Z":fmt(z),"COEFFICIENT_P":fmt(p),"M_TOTAL_M550":fmt(mtotal),"ANNOTATION_M_M550":fmt(ma),"TAU_STAR_CORRECT_M550":fmt(c*sd*mtotal/h2),"TAU_STAR_SE_CORRECT_M550":fmt(se*sd*mtotal/h2),"ILL_CONDITIONED_WARNING":ill,"OVERLAP_DERIVED_ENRICHMENT":"NOT_REPORTED_PHASE12A_REAUDIT_INVALIDATED"})
  q10=bh([float(x["COEFFICIENT_P"]) for x in new])
  for r,q in zip(new,q10):r["BH_FDR_EXPLORATORY_10"]=fmt(q);r["DECISION"]="FDR_SUPPORTED" if q<.05 else ("NOMINAL_ONLY" if float(r["COEFFICIENT_P"])<.05 else "NOT_SUPPORTED")
  write_tsv(OUT/"PHASE17_S_LDSC_RESULTS.tsv",new)
  old=read_tsv(P12B/"PHASE12A_RESULT_RECONCILIATION.tsv")
  old=[x for x in old if x["DISEASE"] in DISEASES]
  combined=[]
  for x in old:combined.append({"FAMILY":"ORIGINAL_FORMAL_30","DISEASE":x["DISEASE"],"CELL_ID":x["CELL_ID"],"COEFFICIENT_P":x["COEFFICIENT_P_RECALC"],"ORIGINAL_BH_FDR_30_PRESERVED":x["BH_FDR_30_RECALC"],"BH_FDR_EXPLORATORY_10":"NA"})
  for x in new:combined.append({"FAMILY":x["FAMILY"],"DISEASE":x["DISEASE"],"CELL_ID":x["CELL_ID"],"COEFFICIENT_P":x["COEFFICIENT_P"],"ORIGINAL_BH_FDR_30_PRESERVED":"NA","BH_FDR_EXPLORATORY_10":x["BH_FDR_EXPLORATORY_10"]})
  q40=bh([float(x["COEFFICIENT_P"]) for x in combined])
  for x,q in zip(combined,q40):x["BH_FDR_COMBINED_40_SENSITIVITY"]=fmt(q)
  write_tsv(OUT/"PHASE17_BH_CORRECTION.tsv",combined)
  counts=Counter(x["retained_factors"] for x in models);supported=FORMAL+CELLS
  audit=f"""# Phase17 read-only audit report\n\n## Frozen inventories\n\nThe Phase9 Swedish all-cell ledger contains 124 models: {counts['0']} retained zero factors, {counts['1']} retained one, {counts['2']} retained two, and {counts['3']} retained three. Factor retention is a model property, not a positive biological result. The complete row-level inventory is `PHASE17_CELL_INVENTORY.tsv`.\n\nThe formal HUNT family contains 18 cells: {', '.join(sorted(formal))}. The four post-Stage-A exploratory cells are {', '.join(exp)}. Exactly eight cells have formal or exploratory cell-level cross-cohort support: {', '.join(supported)}.\n\nPhase12A disease S-LDSC used exactly the six formal SUPPORTED annotations: {', '.join(FORMAL)}. ET and ITL6_2_2 were not included in Phase12A disease S-LDSC. Phase12B generated cell-local factor-score inputs for both but explicitly paused its 124-cell disease extension; this was not a disease result. Phase13-16 did not add ET/ITL6_2_2 disease S-LDSC or molecular-QTL tests. Existing six-cell positional overlaps and Phase15/16 eQTL work used different annotation/cell systems and do not constitute these tests.\n\n## Eligibility decision\n\nBoth target annotations recovered the exact frozen CATlas BED hashes and GRCh38 source-native cCRE definition without a new threshold. They passed 22-chromosome annotation/LD-score construction, HM3 overlap, baselineLD keyed-row coverage, and regression completion gates. Therefore the only eligible extension under the unchanged design is the prespecified 2 cells x 5 diseases marginal S-LDSC family. Phase12B GLOBAL Top5, conjFDR, MR and molecular-QTL analyses are not extended.\n"""
  (OUT/"PHASE17_AUDIT_REPORT.md").write_text(audit,encoding="utf-8")
  methods="""# Phase17 methods\n\nET and ITL6_2_2 were prespecified from the frozen Phase9 post-Stage-A exploratory HUNT family. Their unchanged CATlas 2023 GRCh38 cell-type-assigned accessible cCRE BED files were converted to binary annotations with the same Phase12A builder. Custom LD scores used LDSC v1.0.1, GRCh38 1000 Genomes EUR PLINK files, a 1-cM window, and HM3/no-MHC regression SNPs. Each model included baselineLD v2.2 and one custom annotation, used the original Phase12A munged SCZ, BD, AD, CD or UC summary statistics, the same weights and frequency files, and `--overlap-annot --print-coefficients`.\n\nInference uses the custom-annotation coefficient, its standard error and two-sided normal P value. Following the Phase12B re-audit, overlap-derived Prop_SNPs, Prop_h2 and enrichment are not reported. Tau-star uses `tau * sqrt[p(1-p)] * M / h2`, where p is annotation M_5_50 / total M_5_50 and M=5,956,110. BH is applied separately across the ten new tests. The frozen original 30-test q values are preserved unchanged; a newly calculated 40-test BH is sensitivity only. Nominal P<0.05 without exploratory-family q<0.05 is NOMINAL_ONLY.\n"""
  (OUT/"PHASE17_METHODS.md").write_text(methods,encoding="utf-8")
  hits=[x for x in new if x["DECISION"]=="FDR_SUPPORTED"];nom=[x for x in new if x["DECISION"]=="NOMINAL_ONLY"]
  decision=f"""# Phase17 final scientific decision\n\nAll ten prespecified models completed and both annotations were technically eligible. Exploratory-family FDR-supported results: {', '.join(x['DISEASE']+'–'+x['CELL_LABEL'] for x in hits) if hits else 'none'}. Nominal-only results: {', '.join(x['DISEASE']+'–'+x['CELL_LABEL'] for x in nom) if nom else 'none'}. Every other completed association is NOT_SUPPORTED; none is NOT_ELIGIBLE. Exact coefficients, P values, q values, diagnostics and corrected tau-star are in the result tables.\n\nThese exploratory annotations {('add '+str(len(hits))+' FDR-supported disease-related regulatory contexts beyond the original six-cell family') if hits else 'do not add an FDR-supported disease-related context beyond the original six-cell family'}. The combined 40-test BH is only a sensitivity analysis and does not revise the original Phase12A family. Results remain annotation-level conditional heritability associations; they do not establish a causal cell type, replicated SNP effects, a microbiome-mediated mechanism, or disease relevance of retained GFA factors by themselves.\n"""
  (OUT/"PHASE17_FINAL_SCIENTIFIC_DECISION.md").write_text(decision,encoding="utf-8")
  inputs=[]
  for p in [P9/"01_CELL_ATLAS_SOURCE_MANIFEST.tsv",P9/"04_ALL_CELL_GFA_MODEL_SUMMARY.tsv",P9/"STAGEB_HUNT_VALIDATION_MANIFEST.tsv",P9/"STAGEB_POST_STAGEA_EXPLORATORY_RESULTS.tsv",P12B/"PHASE12A_RESULT_RECONCILIATION.tsv"]+[Path(meta[c]["bed_path"]) for c in CELLS]+[P12/"munged"/(d+".sumstats.gz") for d in DISEASES]:
    inputs.append({"path":str(p),"bytes":p.stat().st_size,"sha256":sha(p)})
  write_tsv(OUT/"PHASE17_INPUT_SHA256_MANIFEST.tsv",inputs)
  outputs=[]
  for p in sorted(OUT.rglob("*")):
    if p.is_file() and p.name not in {"PHASE17_OUTPUT_SHA256_MANIFEST.tsv"}:outputs.append({"path":str(p),"bytes":p.stat().st_size,"sha256":sha(p)})
  write_tsv(OUT/"PHASE17_OUTPUT_SHA256_MANIFEST.tsv",outputs)
if __name__=="__main__":main()
