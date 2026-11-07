#!/usr/bin/env python3
import os
from pathlib import Path
import csv, hashlib, math

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT=ROOT/"PHASE18_MANUSCRIPT_EVIDENCE_FINAL"
P9=ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"
P12=ROOT/"PHASE12A_CELL_DISEASE_SLD SC"
P12B=ROOT/"PHASE12B_MICROBIOME_INFORMED_DISEASE_INTEGRATION"
P14=ROOT/"PHASE14_TARGETED_VALIDATION"/"B_SCZ_CD_UC_LOCI"
P15=ROOT/"PHASE15_CELL_MOLECULAR_AND_MANUSCRIPT"
P16=ROOT/"PHASE16_SCZ_THREE_WAY_SIGNAL_AUDIT"
P17=ROOT/"PHASE17_EXPLORATORY_CELL_DISEASE_EXTENSION"
MW=ROOT/"MANUSCRIPT_WRITING_HANDOFF"

def rd(p):
  with open(p,encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def wr(p,rows,fields=None):
  if not rows: raise ValueError("no rows for "+str(p))
  fields=fields or list(rows[0])
  with open(p,"w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore",lineterminator="\n");w.writeheader();w.writerows(rows)
def bh(vals):
  order=sorted(range(len(vals)),key=lambda i:vals[i]);q=[1.0]*len(vals);run=1.0;n=len(vals)
  for rank in range(n,0,-1):
    i=order[rank-1];run=min(run,vals[i]*n/rank);q[i]=run
  return q
def fmt(x):return "%.12g"%float(x)
def sha(p):
  h=hashlib.sha256()
  with open(p,"rb") as f:
    for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
  return h.hexdigest()

inputs=[
 P9/"04_ALL_CELL_GFA_MODEL_SUMMARY.tsv",P9/"05_CELL_GFA_GLOBAL_ALIGNMENT.tsv",P9/"STAGEB_HUNT_VALIDATION_MANIFEST.tsv",
 P9/"STAGEB_18_CELL_VALIDATION_RESULTS.tsv",P9/"06_CELL_GFA_CROSSCOHORT_COMPOSITION.tsv",P9/"STAGEB_POST_STAGEA_EXPLORATORY_RESULTS.tsv",
 P9/"STAGEB_POST_STAGEA_EXPLORATORY_MANIFEST.tsv",P9/"STAGEB_HUNT_FACTOR_MATCHING_AUDIT.tsv",P9/"STAGEB_HUNT_RBO_NULL_CALIBRATION.tsv",
 P9/"PHASE9C_CONSERVATIVE_NULL_RESULTS.tsv",P9/"PHASE9C_GLOBAL_BASELINE_COMPARISON.tsv",P9/"PHASE9C_PREFERENCE_MAGNITUDE_AUDIT.tsv",
 P9/"PHASE9C_FINAL_SCIENTIFIC_DECISION.md",P9/"PHASE9_CORE_RESULTS_DECISION_PRE_STORYLINE.md",P9/"PHASE9_STAGEB_EXECUTION_FREEZE.md",
 P9/"PHASE9_REVISED_EXECUTION_MANIFEST.tsv",P9/"scripts"/"21_stageb_crosscohort_validation.R",
 P12/"PHASE12A_FINAL_SCIENTIFIC_REPORT.md",P12/"PHASE12A_36_PRIMARY_RESULTS.tsv",P12B/"PHASE12A_RESULT_RECONCILIATION.tsv",
 P17/"PHASE17_FINAL_SCIENTIFIC_DECISION.md",P17/"PHASE17_S_LDSC_RESULTS.tsv",P17/"PHASE17_BH_CORRECTION.tsv",P17/"PHASE17_ANNOTATION_QC.tsv",P17/"PHASE17_METHODS.md",
 P14/"PHASE14B_INTEGRATED_EVIDENCE_MATRIX.tsv",P14/"PHASE14B_CORRECTED_SPECIES_LD_RESULTS.tsv",P14/"PHASE14B_COLOCALIZATION_RESULTS.tsv",P14/"PHASE14B_HUNT_LOOKUP.tsv",
 P15/"SCZ_CHR3_BRAIN_QTL_RESULTS.tsv",P15/"CD_FUT2_QTL_REFERENCE.tsv",P15/"UC_CHR5_GUT_QTL_RESULTS.tsv",P15/"CELL_MOLECULAR_INTEGRATION_REPORT.md",
 P16/"PAIRWISE_COLOC_REANALYSIS.tsv",P16/"ALLELE_DIRECTION_CONCORDANCE.tsv",P16/"THREE_WAY_COLOC_ELIGIBILITY.tsv",P16/"PHASE16_FINAL_SCIENTIFIC_REPORT.md",
 MW/"RESULTS_EVIDENCE_LEDGER.tsv",MW/"FIGURE_AND_SUPPLEMENTARY_TABLE_INVENTORY.md"]

def hunt():
  models=rd(P9/"04_ALL_CELL_GFA_MODEL_SUMMARY.tsv")
  assert len(models)==124
  counts={k:sum(x["retained_factors"]==k for x in models) for k in ["0","1","2","3"]}
  dirs=rd(P9/"06_CELL_GFA_CROSSCOHORT_COMPOSITION.tsv");assert len(dirs)==44
  q44=bh([float(x["empirical_p"]) for x in dirs])
  out=[]
  for x,q in zip(dirs,q44):
    y=dict(x);y["BH_q_UNIFIED_44"]=fmt(q);y["composition_validated_UNIFIED_44"]=str(x["direction_evaluable"]=="TRUE" and q<=.10).upper()
    y["direction_supported_UNIFIED_44"]=str(x["direction_evaluable"]=="TRUE" and q<=.10 and x["directionally_consistent"]=="TRUE").upper()
    y["historical_family_preserved"]=x["BH_family"];y["historical_BH_q_preserved"]=x["BH_q"]
    out.append(y)
  wr(OUT/"PHASE18_HUNT_44_TEST_FINAL.tsv",out)
  formal={x["annotation_id"] for x in rd(P9/"STAGEB_HUNT_VALIDATION_MANIFEST.tsv")}
  explore={x["annotation_id"]:x for x in rd(P9/"STAGEB_POST_STAGEA_EXPLORATORY_MANIFEST.tsv")}
  historical={x["annotation_id"]:x["validation_status"] for x in rd(P9/"STAGEB_18_CELL_VALIDATION_RESULTS.tsv")}
  historical.update({x["annotation_id"]:x["validation_status"] for x in rd(P9/"STAGEB_POST_STAGEA_EXPLORATORY_RESULTS.tsv")})
  status=[]
  for cell in sorted(set(x["annotation_id"] for x in out)):
    rr=[x for x in out if x["annotation_id"]==cell];both=all(x["direction_evaluable"]=="TRUE" for x in rr)
    pref=rr[0]["preference_concordant"]=="TRUE";supported=[x for x in rr if x["direction_supported_UNIFIED_44"]=="TRUE"]
    if both and pref and supported:s="COMPLETE_SUPPORT"
    elif (not both) and supported:s="PARTIAL_SINGLE_DIRECTION_EVIDENCE"
    elif both:s="NOT_SUPPORTED"
    else:s="NOT_EVALUABLE"
    source="ORIGINAL_FORMAL_18" if cell in formal else "POST_STAGE_A_EXPLORATORY_4"
    fr=next(x for x in models if x["annotation_id"]==cell)
    er=explore.get(cell,{})
    status.append({"annotation_id":cell,"organ":rr[0]["organ"],"panel_origin":source,"swedish_retained_factors":fr["retained_factors"],"F1_evaluable":next(x for x in rr if x["direction"]=="F1")["direction_evaluable"],"F2_evaluable":next(x for x in rr if x["direction"]=="F2")["direction_evaluable"],"preference_concordant":rr[0]["preference_concordant"],"supported_directions":";".join(x["direction"] for x in supported) or "NONE","unified_status":s,"historical_status":historical[cell],"selection_frozen_at":er.get("frozen_at",next((x["selection_rules_frozen_at"] for x in rd(P9/"STAGEB_HUNT_VALIDATION_MANIFEST.tsv") if x["annotation_id"]==cell),"")),"pre_HUNT_result_freeze_evidence":("Manifest frozen 2026-10-30 10:14 with HUNT status deferred/0 attempts; result tables generated 2026-10-31 08:01" if cell in explore else "Formal rules frozen 2026-10-24; StageB panel frozen before new HUNT execution"),"disease_S_LDSC_available":str(cell in {"brain__BNGA","brain__SNCG_1","brain__SMC","gut__Macrophage","gut__Early_Colonocyte","gut__NaiveB","brain__ET","brain__ITL6_2_2"}).upper()})
  wr(OUT/"PHASE18_UNIFIED_22_CELL_STATUS.tsv",status)
  return counts,status,out

def disease(status):
  supported={x["annotation_id"] for x in status if x["unified_status"]=="COMPLETE_SUPPORT"}
  old=rd(P12B/"PHASE12A_RESULT_RECONCILIATION.tsv");old=[x for x in old if x["DISEASE"] in {"SCZ","BD","AD","CD","UC"}]
  assert len(old)==30
  new=rd(P17/"PHASE17_S_LDSC_RESULTS.tsv");assert len(new)==10
  rows=[]
  for x in old:
    p=float(x["COEFFICIENT_P_RECALC"]);c=float(x["COEFFICIENT"])
    rows.append({"analysis_family":"HISTORICAL_FORMAL_30","cell_annotation":x["CELL_ID"],"disease":x["DISEASE"],"coefficient":x["COEFFICIENT"],"SE":x["COEFFICIENT_SE"],"P":fmt(p),"historical_BH_q":x["BH_FDR_30_RECALC"],"unified_BH_q_40":"","direction":"POSITIVE" if c>0 else "NEGATIVE","QC_status":"PASS_COEFFICIENT_REAUDIT","tau_star_M550":x["TAU_STAR_CORRECT_M550"],"unified_HUNT_status":"COMPLETE_SUPPORT" if x["CELL_ID"] in supported else "NOT_COMPLETE_SUPPORT","enrichment_status":"INVALID_NOT_USED"})
  for x in new:
    p=float(x["COEFFICIENT_P"]);c=float(x["COEFFICIENT"])
    rows.append({"analysis_family":"HISTORICAL_EXPLORATORY_10","cell_annotation":x["CELL_ID"],"disease":x["DISEASE"],"coefficient":x["COEFFICIENT"],"SE":x["COEFFICIENT_SE"],"P":fmt(p),"historical_BH_q":x["BH_FDR_EXPLORATORY_10"],"unified_BH_q_40":"","direction":"POSITIVE" if c>0 else "NEGATIVE","QC_status":"PASS" if x["ILL_CONDITIONED_WARNING"]=="NO" else "PASS_WITH_WARNING","tau_star_M550":x["TAU_STAR_CORRECT_M550"],"unified_HUNT_status":"COMPLETE_SUPPORT" if x["CELL_ID"] in supported else "NOT_COMPLETE_SUPPORT","enrichment_status":"INVALID_NOT_USED"})
  q=bh([float(x["P"]) for x in rows])
  for x,v in zip(rows,q):x["unified_BH_q_40"]=fmt(v);x["unified_decision_40"]="FDR_SUPPORTED" if v<.05 else ("NOMINAL_ONLY" if float(x["P"])<.05 else "NOT_SUPPORTED")
  phase17_q={(x["DISEASE"],x["CELL_ID"]):float(x["BH_FDR_COMBINED_40_SENSITIVITY"]) for x in rd(P17/"PHASE17_BH_CORRECTION.tsv")}
  assert len(phase17_q)==40
  assert all(abs(float(x["unified_BH_q_40"])-phase17_q[(x["disease"],x["cell_annotation"])])<1e-9 for x in rows)
  wr(OUT/"PHASE18_DISEASE_FINAL_RESULTS.tsv",rows)
  bhr=[{"cell_annotation":x["cell_annotation"],"disease":x["disease"],"P":x["P"],"historical_family":x["analysis_family"],"historical_BH_q_preserved":x["historical_BH_q"],"unified_BH_q_40":x["unified_BH_q_40"],"unified_decision_40":x["unified_decision_40"],"independent_verification":"MATCH_PHASE17_BH_TABLE"} for x in rows]
  wr(OUT/"PHASE18_DISEASE_UNIFIED_BH.tsv",bhr)
  covered={x["cell_annotation"] for x in rows};cov=[]
  for x in status:
    if x["unified_status"]=="COMPLETE_SUPPORT":
      miss=sorted(set(["SCZ","BD","AD","CD","UC"])-{r["disease"] for r in rows if r["cell_annotation"]==x["annotation_id"]})
      cov.append({"cell_annotation":x["annotation_id"],"unified_HUNT_status":x["unified_status"],"disease_models_expected":5,"disease_models_available":5-len(miss),"missing_diseases":";".join(miss) or "NONE","additional_workload_models":len(miss),"coverage_decision":"COMPLETE" if not miss else "MISSING_DO_NOT_RUN_AUTOMATICALLY"})
  wr(OUT/"PHASE18_DISEASE_COVERAGE_AUDIT.tsv",cov)
  assert len(supported)==8 and len(cov)==8 and sum(int(x["additional_workload_models"]) for x in cov)==0
  return rows,cov

def locus():
  corr=rd(P14/"PHASE14B_CORRECTED_SPECIES_LD_RESULTS.tsv")
  coloc=rd(P14/"PHASE14B_COLOCALIZATION_RESULTS.tsv")
  hunt=rd(P14/"PHASE14B_HUNT_LOOKUP.tsv")
  p16=rd(P16/"PAIRWISE_COLOC_REANALYSIS.tsv")
  alleles=rd(P16/"ALLELE_DIRECTION_CONCORDANCE.tsv")
  elig=rd(P16/"THREE_WAY_COLOC_ELIGIBILITY.tsv")
  fut2=rd(P15/"CD_FUT2_QTL_REFERENCE.tsv");uc=rd(P15/"UC_CHR5_GUT_QTL_RESULTS.tsv")
  def cr(region,species):return next(x for x in corr if x["REGION_ID"]==region and x["SPECIES"]==species)
  def co(region,species):return next(x for x in coloc if x["REGION_ID"]==region and x["SPECIES"]==species)
  def hu(region,species):return next(x for x in hunt if x["REGION_ID"]==region and x["SPECIES"]==species)
  rows=[]
  def add(locus,item,evidence,value,prior="NA",universe="NA",placement="SUPPLEMENTARY",claim="",caution="",source=""):
    rows.append({"locus":locus,"evidence_item":item,"evidence_type":evidence,"value":value,"default_prior":prior,"SNP_universe":universe,"manuscript_placement":placement,"supported_claim":claim,"prohibited_overclaim":caution,"source_file":source})
  scz=cr("SCZ_3_17582158_18082158","Pullichristensenella excrementipullorum")
  add("SCZ_chr3","disease_and_species_leads","regional+LD",f"disease rs2337387 G/A beta={scz['DISEASE_LEAD_BETA']} P={scz['DISEASE_LEAD_P']}; species lead {scz['SPECIES_STRONGEST_KEY']} beta={scz['SPECIES_STRONGEST_BETA']} P={scz['SPECIES_STRONGEST_P']}; r2={scz['PAIRED_LEAD_TO_SPECIES_STRONGEST_R2']}",placement="MAIN_TEXT",claim="High-LD regional SCZ/species association in Swedish data",caution="Not independent HUNT replication or mediation",source="PHASE14B_CORRECTED_SPECIES_LD_RESULTS.tsv")
  c=co("SCZ_3_17582158_18082158","Pullichristensenella excrementipullorum")
  add("SCZ_chr3","SCZ_species_coloc","coloc.abf",f"PP3={c['PP_H3']}; PP4={c['PP_H4']}; top={c['COLOC_TOP_SNP']} SNP.PP.H4={c['COLOC_TOP_SNP_PP_H4']}","p1=1e-4;p2=1e-4;p12=1e-5",c["COMMON_SNP_N"],"MAIN_TEXT","Default-prior pairwise shared-signal support","Not proof of a causal variant", "PHASE14B_COLOCALIZATION_RESULTS.tsv")
  add("SCZ_chr3","SCZ_species_prior_sensitivity","coloc.abf prior sensitivity",f"PP4(p12=1e-6)={c['PP_H4_P12_1E6']}; PP4(p12=1e-4)={c['PP_H4_P12_1E4']}","p1=p2=1e-4;alternative p12",c["COMMON_SNP_N"],"SUPPLEMENTARY","Prior dependence is explicit","Do not promote sensitivity maximum to the main estimate","PHASE14B_COLOCALIZATION_RESULTS.tsv")
  for cell in ["Excitatory.neurons","OPCs...COPs"]:
    for comp in ["SCZ_vs_TBC1D5_eQTL","species_vs_TBC1D5_eQTL"]:
      z=next(x for x in p16 if x["cell_type"]==cell and x["comparison"]==comp and x["analysis_scope"]=="three_way_common_501" and x["p12"]=="1e-05")
      add("SCZ_chr3",cell+"_"+comp,"pairwise coloc common universe",f"PP3={z['PP_H3']}; PP4={z['PP_H4']}; top={z['top_snp']}; rs2337387 SNP.PP.H4={z['rs2337387_SNP_PP_H4']}","p1=1e-4;p2=1e-4;p12=1e-5","501 common non-palindromic SNPs","MAIN_TEXT" if cell=="Excitatory.neurons" else "SUPPLEMENTARY","Pairwise molecular-QTL context in a harmonized SNP universe","Not formal three-way colocalization", "PAIRWISE_COLOC_REANALYSIS.tsv")
      sens=[a for a in p16 if a["cell_type"]==cell and a["comparison"]==comp and a["analysis_scope"]=="three_way_common_501" and a["p12"] in {"1e-06","1e-04"}]
      add("SCZ_chr3",cell+"_"+comp+"_prior_sensitivity","pairwise coloc prior sensitivity","; ".join("p12="+a["p12"]+":PP4="+a["PP_H4"] for a in sens),"p1=p2=1e-4;alternative p12","501 common non-palindromic SNPs","SUPPLEMENTARY","Prior dependence is explicit","Do not select priors to manufacture concordance","PAIRWISE_COLOC_REANALYSIS.tsv")
  for cell in ["Excitatory.neurons","OPCs...COPs"]:
    z=next(x for x in alleles if x["cell_type"]==cell and x["variant_roles"]=="prespecified_scz_lead")
    add("SCZ_chr3",cell+"_rs2337387_direction","allele direction",f"SCZ risk allele {z['scz_risk_effect_allele']}: SCZ beta={z['scz_beta_per_risk_allele']}; species beta={z['species_beta_per_scz_risk_allele']}; TBC1D5 eQTL beta={z['tbc1d5_eqtl_beta_per_scz_risk_allele']}",universe="rs2337387 harmonized allele",placement="MAIN_TEXT" if cell=="Excitatory.neurons" else "SUPPLEMENTARY",claim="Cell-context-dependent TBC1D5 direction",caution="Direction does not establish mediation",source="ALLELE_DIRECTION_CONCORDANCE.tsv")
  h=hu("SCZ_3_17582158_18082158","Pullichristensenella excrementipullorum")
  add("SCZ_chr3","HUNT_species_lookup","descriptive external lookup",f"lead P={h['HUNT_LEAD_P']}; direction_concordant={h['LEAD_DIRECTION_CONCORDANT']}",placement="SUPPLEMENTARY",claim="No HUNT support at the disease lead",caution="Do not call replicated",source="PHASE14B_HUNT_LOOKUP.tsv")
  add("SCZ_chr3","formal_three_way_eligibility","HyPrColoc/moloc eligibility","NOT_RUN_INPUT_ASSUMPTIONS_NOT_VERIFIABLE",universe="501 common SNPs",placement="MAIN_TEXT",claim="Pairwise evidence only",caution="No three-way causal colocalization",source="THREE_WAY_COLOC_ELIGIBILITY.tsv")
  cd=cr("CD_19_48453417_48953417","Acutalibacter ornithocaccae");cc=co("CD_19_48453417_48953417","Acutalibacter ornithocaccae");fh=fut2[0]
  add("CD_FUT2","disease_species_LD","regional+LD",f"rs601338 A/G beta={cd['DISEASE_LEAD_BETA']} P={cd['DISEASE_LEAD_P']}; species lead {cd['SPECIES_STRONGEST_KEY']} P={cd['SPECIES_STRONGEST_P']}; r2={cd['PAIRED_LEAD_TO_SPECIES_STRONGEST_R2']}",placement="MAIN_TEXT",claim="Known FUT2 host-microbiome reference",caution="Not novel",source="PHASE14B_CORRECTED_SPECIES_LD_RESULTS.tsv")
  add("CD_FUT2","CD_species_coloc","coloc.abf",f"PP4={cc['PP_H4']}; PP3={cc['PP_H3']}","p1=1e-4;p2=1e-4;p12=1e-5",cc["COMMON_SNP_N"],"MAIN_TEXT","Pairwise sensitivity support","Original species MAF/INFO/SNP-N absent", "PHASE14B_COLOCALIZATION_RESULTS.tsv")
  add("CD_FUT2","CD_species_prior_sensitivity","coloc.abf prior sensitivity",f"PP4(p12=1e-6)={cc['PP_H4_P12_1E6']}; PP4(p12=1e-4)={cc['PP_H4_P12_1E4']}","p1=p2=1e-4;alternative p12",cc["COMMON_SNP_N"],"SUPPLEMENTARY","Prior dependence is explicit","Known reference does not remove input limitations","PHASE14B_COLOCALIZATION_RESULTS.tsv")
  add("CD_FUT2","CEDAR2_FUT2","disease/species eQTL context",f"node={fh['cell_tree_id']}; rs601338 effect allele A; eQTL beta={fh['rs601338_eqtl_beta']} P={fh['rs601338_eqtl_p']}; CD-eQTL PP4={fh['cd_eqtl_pp4_default']}; species-eQTL PP4={fh['best_species_eqtl_pp4_default']}","p1=1e-4;p2=1e-4;p12=1e-5",fh["harmonized_coloc_snp_n"],"MAIN_TEXT","Strong known FUT2 regulatory reference","Not a novel microbiome-mediated mechanism", "CD_FUT2_QTL_REFERENCE.tsv")
  h=hu("CD_19_48453417_48953417","Acutalibacter ornithocaccae");add("CD_FUT2","HUNT_species_lookup","descriptive external lookup",f"lead P={h['HUNT_LEAD_P']}; direction_concordant={h['LEAD_DIRECTION_CONCORDANT']}",placement="SUPPLEMENTARY",claim="Concordant descriptive external evidence",caution="Not a preregistered locus replication family",source="PHASE14B_HUNT_LOOKUP.tsv")
  for sp in ["Enterocloster excrementigallinarum","NSJ-51 sp900549795"]:
    z=cr("UC_5_308192_808192",sp);cz=co("UC_5_308192_808192",sp)
    add("UC_chr5",sp,"regional LD + coloc",f"species lead={z['SPECIES_STRONGEST_KEY']} P={z['SPECIES_STRONGEST_P']} r2={z['PAIRED_LEAD_TO_SPECIES_STRONGEST_R2']}; PP4={cz['PP_H4']}","p1=1e-4;p2=1e-4;p12=1e-5",cz["COMMON_SNP_N"],"SUPPLEMENTARY","NSJ-51, not Enterocloster, has high-LD/default-prior coloc support","Do not transfer evidence between species", "PHASE14B_CORRECTED_SPECIES_LD_RESULTS.tsv;PHASE14B_COLOCALIZATION_RESULTS.tsv")
    add("UC_chr5",sp+"_prior_sensitivity","coloc.abf prior sensitivity",f"PP4(p12=1e-6)={cz['PP_H4_P12_1E6']}; PP4(p12=1e-4)={cz['PP_H4_P12_1E4']}","p1=p2=1e-4;alternative p12",cz["COMMON_SNP_N"],"SUPPLEMENTARY","Prior dependence is explicit","Do not substitute sensitivity PP4 for default-prior evidence","PHASE14B_COLOCALIZATION_RESULTS.tsv")
  for z in uc:add("UC_chr5","EXOC3_"+z["cell_tree_id"],"CEDAR-2 local eQTL",f"lead={z['permutational_lead_variant']} nominal P={z['permutational_nominal_p']} q={z['permutational_q']}; conditional region SNPs={z['conditional_rank0_snp_n_in_region']}; harmonized SNPs={z['harmonized_disease_eqtl_snp_n']}; status={z['coloc_abf_status']}",universe=z["harmonized_disease_eqtl_snp_n"],placement="SUPPLEMENTARY",claim="Local EXOC3 regulatory context only",caution="Do not claim EXOC3 colocalization",source="UC_CHR5_GUT_QTL_RESULTS.tsv")
  wr(OUT/"PHASE18_LOCUS_MOLECULAR_EVIDENCE.tsv",rows)
  return rows

def docs(counts,status,disease_rows,cov,locus_rows):
  complete=[x["annotation_id"] for x in status if x["unified_status"]=="COMPLETE_SUPPORT"]
  partial=[x["annotation_id"] for x in status if x["unified_status"]=="PARTIAL_SINGLE_DIRECTION_EVIDENCE"]
  hits=[x for x in disease_rows if x["unified_decision_40"]=="FDR_SUPPORTED"]
  ev=f"""# Phase18 Results 4–6 evidence map\n\n## Results 4 — cell-regulatory discovery and HUNT evaluation\n\n**Result.** Swedish models completed for 124 annotations, retaining 0/1/2/3 factors in {counts['0']}/{counts['1']}/{counts['2']}/{counts['3']} annotations. Factor retention alone is not biological support. The original 18-cell HUNT panel and four post-Stage-A Swedish-derived candidates remain distinct selection strata. The latter manifest was frozen on 2026-10-30 at 10:14 with HUNT status deferred and zero attempts; unified result files were generated on 2026-10-31 at 08:01. Thus local records support pre-HUNT-result finalization, although no external version-control commit proves authorship chronology beyond these immutable manifests/timestamps.\n\nAcross the unified 44 planned directions, non-evaluable tests retained P=1 and BH used all 44. Complete support was retained for eight cells: {', '.join(complete)}. Partial single-direction evidence occurred for {', '.join(partial)}. These are not eight equivalently selected replicated cells: six came from the original formal panel and ET/ITL6_2_2 from the separate exploratory panel.\n\n**Placement.** Main text Results 4; Figure 3 cell-regulatory panel; Phase18 Tables `PHASE18_HUNT_44_TEST_FINAL.tsv` and `PHASE18_UNIFIED_22_CELL_STATUS.tsv`; historical workbook S15–S25 for provenance. Background comparisons remain sensitivity only because all 124 matched backgrounds failed full functional-covariate balance.\n\n## Results 5 — unified disease S-LDSC\n\nAll eight completely supported cells have five existing disease models; additional workload is 0. The manuscript-integrated family should therefore be the 40 coefficient tests across eight cells and five diseases. Historical 30-test and 10-test q values remain traceability columns, but final inferential wording should use unified q40. Unified FDR-supported coefficients are: {', '.join(x['disease']+'–'+x['cell_annotation']+' ('+x['direction'].lower()+', q='+x['unified_BH_q_40']+')' for x in hits)}. Positive and negative coefficients are conditional baselineLD-adjusted directions, not marginal enrichment or depletion.\n\n**Placement.** Main text Results 5; rebuilt Figure 3 coefficient/tau-star heatmap; `PHASE18_DISEASE_FINAL_RESULTS.tsv`, `PHASE18_DISEASE_UNIFIED_BH.tsv`, and coverage audit. Phase12A marginal-enrichment Figure 2 is DO_NOT_USE.\n\n## Results 6 — targeted loci and molecular QTL context\n\n**SCZ chr3.** Swedish SCZ–Pullichristensenella evidence has high lead-proxy LD and default-prior pairwise coloc support; TBC1D5 disease/eQTL and actual species/eQTL comparisons provide partial pairwise molecular context. Excitatory-neuron and OPC/COP eQTL effects differ in direction. HUNT does not support the species association at rs2337387. HyPrColoc/moloc were ineligible; no three-way colocalization or microbial mediation claim is allowed. Main text with Figure 4.\n\n**CD FUT2.** rs601338, Acutalibacter, pairwise coloc, and CEDAR-2 FUT2 eQTL evidence form a strong known positive reference. The rs601338 A allele decreases FUT2 expression in the leading CEDAR-2 node. This is not a novel discovery. Main text with Figure 5.\n\n**UC chr5.** Enterocloster has a strong regional association but weak lead LD and low default-prior PP4; the high-LD/default-prior coloc evidence belongs to NSJ-51. EXOC3 has local eQTL evidence but only one or two conditional/harmonized SNPs, so coloc is ineligible. Supplementary-only, with the UC contrast panel of Figure 5 if retained.\n\nExact priors, SNP universes, allele directions, source rows and prohibited claims are in `PHASE18_LOCUS_MOLECULAR_EVIDENCE.tsv`.\n"""
  (OUT/"PHASE18_RESULTS_4_6_EVIDENCE_MAP.md").write_text(ev,encoding="utf-8")
  methods="""# Phase18 methods evidence map\n\n## Results 4 methods\n\nSource files: Phase9 all-cell model summary, global alignment, original 18-cell manifest, four-cell exploratory manifest, factor-matching audit, empirical RBO null table and cross-cohort direction table. Recalculation changes only multiplicity bookkeeping: original empirical permutation P values were subjected to BH across all 44 planned directions, with non-evaluable directions fixed at P=1. No permutations or GFA models were rerun. Cell-level complete support required two directions evaluable in both cohorts, preference concordance, and at least one direction with q44<=0.10 plus positive signed cosine and common-Top20 direction concordance>=0.50. A supported evaluable direction in an otherwise two-direction-incomplete cell was labelled partial single-direction evidence.\n\nThe original 18 and four post-Stage-A candidates must remain separate selection strata in Methods. File timestamps and frozen manifest fields support that the four candidates were fixed before their HUNT model attempts/results; this is a local provenance audit, not an external preregistration claim.\n\n## Results 5 methods\n\nThe 40 already completed single-annotation S-LDSC models were harmonized without rerunning LDSC. Coefficients, standard errors and coefficient P values were read from the Phase12B re-audit (original 30) and Phase17 (additional 10). BH across all 40 is the recommended final integrated family because the unified HUNT audit yields the same eight completely supported cells and disease coverage is complete. Historical q30 and q10 are retained. Invalid overlap-derived Prop_SNPs, Prop_h2 and enrichment are excluded; corrected M_5_50 tau-star is descriptive scale information only.\n\n## Results 6 methods\n\nAll locus analyses are read-only harmonizations of Phase14–16. Default pairwise coloc priors are p1=p2=1e-4 and p12=1e-5; alternative p12 values belong only in supplementary sensitivity evidence. Evidence from different SNP universes is shown in separate rows. Disease–eQTL, species–eQTL and disease–species coloc are distinct comparisons. No formal three-trait posterior was calculated.\n"""
  (OUT/"PHASE18_METHODS_EVIDENCE_MAP.md").write_text(methods,encoding="utf-8")
  dec=f"""# Phase18 final scientific decision\n\nThe unified 44-direction BH audit preserves eight completely supported cell annotations: {', '.join(complete)}. This is the same eight-cell set represented by the 40 existing disease S-LDSC models, so no eligible cell-disease model is missing and no additional S-LDSC run is warranted. The final manuscript should use the unified 40-test coefficient family while retaining historical q30/q10 columns.\n\nResults 4 can state cross-cohort support for six formally selected cells plus exploratory support for ET and ITL6_2_2, while reporting partial single-direction and null/non-evaluable results. Results 5 can report unified-FDR coefficient associations, with sign described as conditional direction rather than enrichment/depletion. Results 6 can present SCZ chr3 as an exploratory pairwise concordance candidate, FUT2 as a known positive reference, and UC chr5 as supplementary mixed species-specific evidence.\n\nNo evidence supports calling all 22 cells replicated, claiming complete cell specificity, three-way SCZ colocalization, microbiome mediation, a novel FUT2 mechanism, or EXOC3 colocalization.\n"""
  (OUT/"PHASE18_FINAL_SCIENTIFIC_DECISION.md").write_text(dec,encoding="utf-8")

def manifests():
  wr(OUT/"INPUT_SHA256_MANIFEST.tsv",[{"path":str(p),"bytes":p.stat().st_size,"sha256":sha(p)} for p in inputs])
  files=[]
  for p in sorted(OUT.rglob("*")):
    if p.is_file() and p.name!="OUTPUT_SHA256_MANIFEST.tsv":files.append({"path":str(p),"bytes":p.stat().st_size,"sha256":sha(p)})
  wr(OUT/"OUTPUT_SHA256_MANIFEST.tsv",files)

def main():
  counts,status,h44=hunt();drows,cov=disease(status);lrows=locus();docs(counts,status,drows,cov,lrows);manifests()
if __name__=="__main__":main()
