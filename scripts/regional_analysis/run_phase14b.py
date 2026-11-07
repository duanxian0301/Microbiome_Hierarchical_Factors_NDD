from __future__ import annotations

import csv
import gzip
import itertools
import json
import math
from collections import defaultdict
import os
from pathlib import Path

import numpy as np
from pyliftover import LiftOver


ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
BASE = ROOT / "FULL_SCALE_GFA_FINAL_EVALUATION"
P11 = ROOT / "PHASE11_DISEASE_BIOLOGICAL_INTEGRATION"
P13 = ROOT / "PHASE13_SHARED_STRUCTURE_AND_DISEASE_LOCI" / "B_FOUR_DISEASE_LOCI"
OUT = ROOT / "PHASE14_TARGETED_VALIDATION" / "B_SCZ_CD_UC_LOCI"
INPUTS = OUT / "coloc_inputs"
INPUTS.mkdir(parents=True, exist_ok=True)
REF = ROOT / "GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT" / "CELL_ANNOTATION_C1_EXECUTION" / "reference_download" / "GRCh38" / "plink_files"
RAW = P11 / "02_DISEASE_RESOURCES" / "PHASE11C_RAW"
CHAIN = P11 / "PHASE11E_FUT2_LOCAL_AUDIT" / "hg19ToHg38.over.chain.gz"

REGIONS = [
    dict(region_id="SCZ_3_17582158_18082158", disease="SCZ", chrom=3, start=17582158, end=18082158, lead="rs2337387", lead_bp=17832158, candidate="Pullichristensenella excrementipullorum", gene="BALR6"),
    dict(region_id="SCZ_11_57492822_57992822", disease="SCZ", chrom=11, start=57492822, end=57992822, lead="rs9420", lead_bp=57742822, candidate="Ventrisoma faecale", gene="ENSG00000254732"),
    dict(region_id="UC_5_308192_808192", disease="UC", chrom=5, start=308192, end=808192, lead="rs17497684", lead_bp=558192, candidate="Enterocloster excrementigallinarum", gene="MIR4456"),
    dict(region_id="UC_1_161252368_161752368", disease="UC", chrom=1, start=161252368, end=161752368, lead="rs10800309", lead_bp=161502368, candidate="Dysosmobacter faecalis", gene="FCGR2A"),
    dict(region_id="CD_6_166771837_167271837", disease="CD", chrom=6, start=166771837, end=167271837, lead="rs7749278", lead_bp=167021837, candidate="ER4 sp900317525", gene="ENSG00000285730"),
    dict(region_id="CD_19_48453417_48953417", disease="CD", chrom=19, start=48453417, end=48953417, lead="rs601338", lead_bp=48703417, candidate="Acutalibacter ornithocaccae", gene="FUT2"),
]
REGION_BY_ID = {r["region_id"]: r for r in REGIONS}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def pz(z: float) -> float:
    return math.erfc(abs(float(z))/math.sqrt(2)) if math.isfinite(float(z)) else float("nan")


def canon(chrom: int | str, bp: int, a1: str, a2: str) -> tuple[str, int, str, str]:
    a,b=sorted((a1.upper(),a2.upper()))
    return str(chrom),int(bp),a,b


def parse_key(key: str) -> tuple[str,int,str,str]:
    c,p,a1,a2=key.split(":")
    return canon(c,int(p),a1,a2)


def load_reference():
    keys={}; rs={}; region_rows=defaultdict(list)
    target_chroms={r["chrom"] for r in REGIONS}
    for chrom in sorted(target_chroms):
        with (REF/f"1000G.EUR.hg38.{chrom}.bim").open() as handle:
            for i,line in enumerate(handle):
                c,rid,_,bp,a1,a2=line.split(); bp=int(bp)
                rec=(str(chrom),i,rid,bp,a1.upper(),a2.upper())
                keys[canon(chrom,bp,a1,a2)]=rec; rs[rid]=rec
                for reg in REGIONS:
                    if reg["chrom"]==chrom and reg["start"]<=bp<=reg["end"]:
                        region_rows[reg["region_id"]].append(rec)
    return keys,rs,region_rows


def genotype_reader():
    nind=489; nbyte=(nind+3)//4
    lut=np.empty((256,4),np.float32); code=[0.0,np.nan,1.0,2.0]
    for byte in range(256):
        for j in range(4): lut[byte,j]=code[(byte>>(2*j))&3]
    handles={}
    for chrom in sorted({r["chrom"] for r in REGIONS}):
        h=(REF/f"1000G.EUR.hg38.{chrom}.bed").open("rb")
        if h.read(3)!=b"\x6c\x1b\x01": raise ValueError("BED header")
        handles[str(chrom)]=h
    cache={}
    def geno(rec):
        if rec is None: return None
        key=(rec[0],rec[1])
        if key not in cache:
            h=handles[rec[0]]; h.seek(3+rec[1]*nbyte)
            cache[key]=lut[np.frombuffer(h.read(nbyte),dtype=np.uint8)].reshape(-1)[:nind].copy()
        return cache[key]
    def r2(a,b):
        if a is None or b is None or a[0]!=b[0]: return float("nan")
        x=geno(a); y=geno(b); ok=np.isfinite(x)&np.isfinite(y)
        if ok.sum()<30 or np.std(x[ok])==0 or np.std(y[ok])==0:return float("nan")
        return float(np.corrcoef(x[ok],y[ok])[0,1]**2)
    def maf(rec):
        x=geno(rec)
        if x is None:return float("nan")
        x=x[np.isfinite(x)]
        if not len(x):return float("nan")
        af=float(x.mean()/2); return min(af,1-af)
    return geno,r2,maf,handles


def disease_scz(ref_rs: dict) -> dict[str,list[dict]]:
    target={r["lead"] for r in REGIONS if r["disease"]=="SCZ"}
    for rid,rec in ref_rs.items():
        for reg in REGIONS:
            if reg["disease"]=="SCZ" and rec[0]==str(reg["chrom"]) and reg["start"]<=rec[3]<=reg["end"]: target.add(rid)
    out=defaultdict(list)
    path=RAW/"PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"
    with gzip.open(path,"rt",encoding="latin1") as handle:
        reader=None
        for line in handle:
            if line.startswith("##"): continue
            reader=csv.DictReader(itertools.chain([line],handle),delimiter="\t"); break
        if reader is None: raise ValueError("SCZ header")
        for row in reader:
            if row["ID"] not in target: continue
            rec=ref_rs.get(row["ID"])
            if rec is None: continue
            a1=row["A1"].upper();a2=row["A2"].upper()
            if set((a1,a2))!=set((rec[4],rec[5])):
                ca1=a1.translate(COMPLEMENT);ca2=a2.translate(COMPLEMENT)
                if set((ca1,ca2))==set((rec[4],rec[5])): a1,a2=ca1,ca2
                else: continue
            for reg in REGIONS:
                if reg["disease"]=="SCZ" and rec[0]==str(reg["chrom"]) and reg["start"]<=rec[3]<=reg["end"]:
                    out[reg["region_id"]].append(dict(key=canon(rec[0],rec[3],a1,a2),snp=row["ID"],bp=rec[3],a1=a1,a2=a2,beta=float(row["BETA"]),se=float(row["SE"]),p=float(row["PVAL"]),info=float(row["IMPINFO"]),ncase=float(row["NCAS"]),ncontrol=float(row["NCON"]),source_build="GRCh37 rsID mapped to GRCh38 EUR BIM",maf_source="1000G EUR reference"))
    return out


def disease_ibd(disease: str) -> dict[str,list[dict]]:
    fname={"CD":"cd_build37_40266_20161107.txt.gz","UC":"uc_build37_45975_20161107.txt.gz"}[disease]
    ncase,ncontrol={"CD":(12194,28072),"UC":(12366,33609)}[disease]
    lift=LiftOver(str(CHAIN)); out=defaultdict(list)
    regs=[r for r in REGIONS if r["disease"]==disease]
    with gzip.open(RAW/fname,"rt",encoding="utf-8") as handle:
        for row in csv.DictReader(handle,delimiter="\t"):
            parts=row["MarkerName"].replace("_",":").split(":")
            if len(parts)<4: continue
            try: chrom=int(parts[0]); b37=int(parts[1])
            except ValueError: continue
            possible=[r for r in regs if r["chrom"]==chrom and r["start"]-2_000_000<=b37<=r["end"]+2_000_000]
            if not possible: continue
            a1=row["Allele1"].upper();a2=row["Allele2"].upper()
            if len(a1)!=1 or len(a2)!=1 or a1 not in "ACGT" or a2 not in "ACGT": continue
            hit=lift.convert_coordinate(f"chr{chrom}",b37-1)
            if len(hit)!=1 or hit[0][0]!=f"chr{chrom}": continue
            bp=int(hit[0][1])+1; strand=hit[0][2]
            if strand=="-":a1=a1.translate(COMPLEMENT);a2=a2.translate(COMPLEMENT)
            for reg in possible:
                if reg["start"]<=bp<=reg["end"]:
                    try: beta=float(row["Effect"]);se=float(row["StdErr"]);p=float(row["P.value"])
                    except ValueError: continue
                    out[reg["region_id"]].append(dict(key=canon(chrom,bp,a1,a2),snp=row["MarkerName"],bp=bp,a1=a1,a2=a2,beta=beta,se=se,p=p,info=float("nan"),ncase=ncase,ncontrol=ncontrol,source_build=f"GRCh37 lifted to GRCh38 strand {strand}",maf_source="1000G EUR reference"))
    return out


def build_region_indices(keys: np.ndarray) -> list[np.ndarray]:
    by_chr=defaultdict(list)
    for i,r in enumerate(REGIONS):by_chr[r["chrom"]].append((r["start"],r["end"],i))
    for c in by_chr:by_chr[c].sort()
    out=[[] for _ in REGIONS]; current=None;ptr=0;active=[];intervals=[]
    for idx,raw in enumerate(keys):
        c,p,_,_=str(raw).split(":");c=int(c);p=int(p)
        if c!=current:current=c;ptr=0;active=[];intervals=by_chr.get(c,[])
        while ptr<len(intervals) and intervals[ptr][0]<=p:
            active.append((intervals[ptr][1],intervals[ptr][2]));ptr+=1
        active=[x for x in active if x[0]>=p]
        for _,ri in active:out[ri].append(idx)
    return [np.asarray(x,dtype=np.int64) for x in out]


def load_beds():
    paths={
      "BNGA":BASE/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"annotations"/"brain"/"CATlas2023_BNGA_GRCh38_autosomes.bed.gz",
      "SNCG_1":BASE/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"annotations"/"brain"/"CATlas2023_SNCG_1_GRCh38_autosomes.bed.gz",
      "SMC":BASE/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"annotations"/"brain"/"CATlas2023_SMC_GRCh38_autosomes.bed.gz",
      "Macrophage":BASE/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"annotations"/"gut"/"Zhao2026_Macrophage_GRCh38_autosomes.bed.gz",
      "NaiveB":BASE/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"annotations"/"gut"/"Zhao2026_NaiveB_GRCh38_autosomes.bed.gz",
      "Early_Colonocyte":BASE/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"/"annotations"/"gut"/"Zhao2026_Early_Colonocyte_GRCh38_autosomes.bed.gz",
    }
    target={str(r["chrom"]) for r in REGIONS}; out={}
    for cell,path in paths.items():
        d=defaultdict(list)
        with gzip.open(path,"rt") as h:
            for line in h:
                if not line.strip() or line.startswith("#"):continue
                v=line.rstrip().split("\t");c=v[0].removeprefix("chr")
                if c in target:d[c].append((int(v[1])+1,int(v[2])))
        out[cell]=d
    return out


def in_peak(beds,cell,chrom,bp):
    return any(a<=bp<=b for a,b in beds[cell].get(str(chrom),[]))


def main() -> None:
    manifest=read_tsv(BASE/"inputs_241_manifest.tsv")
    h2=read_tsv(BASE/"BIDIRECTIONAL_SPECIES_LOADING_H2.tsv")
    m_by_acc={r["swedish_accession"]:r for r in manifest}
    species=[]
    for r in h2:
        if r["Swedish_F1_top40"]=="TRUE" or r["Swedish_F2_top40"]=="TRUE":
            species.append(dict(r,primary=(r["Swedish_F1_top20"]=="TRUE" or r["Swedish_F2_top20"]=="TRUE"),sw_path=m_by_acc[r["swedish_accession"]]["swedish_matrix"],hunt_path=m_by_acc[r["swedish_accession"]]["hunt_matrix"]))
    species.sort(key=lambda r:(not r["primary"],int(r["trait_order"])))
    if sum(s["primary"] for s in species)!=31 or len(species)!=63:raise ValueError("species core")

    ref_key,ref_rs,ref_region_rows=load_reference(); geno,r2fun,maffun,handles=genotype_reader()
    disease=defaultdict(list)
    for d in (disease_scz(ref_rs),disease_ibd("CD"),disease_ibd("UC")):
        for k,v in d.items():disease[k].extend(v)
    disease_map={rid:{x["key"]:x for x in rows} for rid,rows in disease.items()}

    first=np.load(manifest[0]["swedish_matrix"],allow_pickle=False); keys=first["snp"]
    region_indices=build_region_indices(keys)
    key_to_idx={}
    for ri,idxs in enumerate(region_indices):
        for idx in idxs:key_to_idx[parse_key(str(keys[idx]))]=int(idx)
    old_candidates={r["REGION_ID"]:r for r in read_tsv(P13/"PHASE13B_PRIORITIZED_CANDIDATES.tsv")}
    rows=[]; hunt_rows=[]; candidate_arrays={}
    candidate_names={r["candidate"] for r in REGIONS}
    for sp in species:
        sw=np.load(sp["sw_path"],allow_pickle=False); hu=np.load(sp["hunt_path"],allow_pickle=False)
        if sp["species"] in candidate_names:
            candidate_arrays[sp["species"]]=(sw["beta"].copy(),sw["se"].copy(),sw["z"].copy())
        for ri,reg in enumerate(REGIONS):
            idxs=region_indices[ri]; zz=sw["z"][idxs]; loc=int(np.nanargmax(np.abs(zz))); idx=int(idxs[loc])
            ckey=parse_key(str(keys[idx])); lead_rec=ref_rs[reg["lead"]]
            lead_key=canon(lead_rec[0],lead_rec[3],lead_rec[4],lead_rec[5]); lead_idx=key_to_idx.get(lead_key)
            lead_beta=float(sw["beta"][lead_idx]) if lead_idx is not None else float("nan")
            lead_se=float(sw["se"][lead_idx]) if lead_idx is not None else float("nan")
            lead_z=float(sw["z"][lead_idx]) if lead_idx is not None else float("nan")
            skey=str(keys[lead_idx]) if lead_idx is not None else "NA"
            dlead=min(disease[reg["region_id"]],key=lambda x:(x["snp"]!=reg["lead"],abs(x["bp"]-reg["lead_bp"])))
            species_a1=skey.split(":")[2] if lead_idx is not None else "NA"; species_a2=skey.split(":")[3] if lead_idx is not None else "NA"
            if dlead["a1"]==species_a1:aligned=lead_beta
            elif dlead["a1"]==species_a2:aligned=-lead_beta
            else:aligned=float("nan")
            paired_r2=r2fun(lead_rec,ref_key.get(ckey))
            old=old_candidates[reg["region_id"]]
            rows.append({"REGION_ID":reg["region_id"],"DISEASE":reg["disease"],"CHR":reg["chrom"],"START":reg["start"],"END":reg["end"],"DISEASE_LEAD_SNP":reg["lead"],"DISEASE_LEAD_BP":reg["lead_bp"],"DISEASE_LEAD_EFFECT_ALLELE":dlead["a1"],"DISEASE_LEAD_OTHER_ALLELE":dlead["a2"],"DISEASE_LEAD_BETA":dlead["beta"],"DISEASE_LEAD_SE":dlead["se"],"DISEASE_LEAD_P":dlead["p"],"SPECIES":sp["species"],"PRIMARY_TOP20_UNION":sp["primary"],"EXPANDED_TOP40_UNION":True,"SPECIES_STRONGEST_KEY":str(keys[idx]),"SPECIES_STRONGEST_BETA":float(sw["beta"][idx]),"SPECIES_STRONGEST_SE":float(sw["se"][idx]),"SPECIES_STRONGEST_Z":float(sw["z"][idx]),"SPECIES_STRONGEST_P":pz(float(sw["z"][idx])),"PAIRED_LEAD_TO_SPECIES_STRONGEST_R2":paired_r2,"LEAD_SPECIES_KEY":skey,"LEAD_SPECIES_ORIGINAL_EFFECT_ALLELE":species_a1,"LEAD_SPECIES_ORIGINAL_BETA":lead_beta,"LEAD_SPECIES_ORIGINAL_SE":lead_se,"LEAD_SPECIES_ORIGINAL_Z":lead_z,"LEAD_SPECIES_ORIGINAL_P":pz(lead_z),"LEAD_SPECIES_BETA_ALIGNED_TO_DISEASE_EFFECT_ALLELE":aligned,"SAME_EFFECT_ALLELE_DIRECTION":"same" if aligned*dlead["beta"]>0 else ("opposite" if aligned*dlead["beta"]<0 else "unresolved"),"PHASE13_BEST_SPECIES":old["BEST_SPECIES"],"PHASE13_REPORTED_MAX_R2":old["MAX_LEAD_TO_SPECIES_CANDIDATE_R2"],"PAIRING_MISMATCH_FOR_PHASE13_BEST_SPECIES":sp["species"]==old["BEST_SPECIES"] and abs(float(old["MAX_LEAD_TO_SPECIES_CANDIDATE_R2"])-paired_r2)>1e-6,"F1_LOADING":sp["Swedish_F1_loading_raw"],"F2_LOADING":sp["Swedish_F2_loading_raw"]})

            hzz=hu["z"][idxs]; hloc=int(np.nanargmax(np.abs(hzz))); hidx=int(idxs[hloc])
            if sp["species"]==reg["candidate"]:
                hlead_z=float(hu["z"][lead_idx]) if lead_idx is not None else float("nan")
                hunt_rows.append({"REGION_ID":reg["region_id"],"DISEASE":reg["disease"],"SPECIES":sp["species"],"LEAD_SNP":reg["lead"],"HUNT_LEAD_BETA":float(hu["beta"][lead_idx]),"HUNT_LEAD_SE":float(hu["se"][lead_idx]),"HUNT_LEAD_Z":hlead_z,"HUNT_LEAD_P":pz(hlead_z),"HUNT_REGION_STRONGEST_KEY":str(keys[hidx]),"HUNT_REGION_STRONGEST_BETA":float(hu["beta"][hidx]),"HUNT_REGION_STRONGEST_SE":float(hu["se"][hidx]),"HUNT_REGION_STRONGEST_Z":float(hu["z"][hidx]),"HUNT_REGION_STRONGEST_P":pz(float(hu["z"][hidx])),"LEAD_TO_HUNT_STRONGEST_R2":r2fun(lead_rec,ref_key.get(parse_key(str(keys[hidx])))),"SWEDISH_LEAD_Z":lead_z,"LEAD_DIRECTION_CONCORDANT":np.sign(lead_z)==np.sign(hlead_z),"INTERPRETATION":"auxiliary lookup; not independent replication unless statistically supported"})
    # Add fixed-region support counts without changing the prespecified threshold.
    for reg in REGIONS:
        q=[x for x in rows if x["REGION_ID"]==reg["region_id"] and x["PRIMARY_TOP20_UNION"]]
        n=sum(float(x["SPECIES_STRONGEST_P"])<1e-4 and float(x["PAIRED_LEAD_TO_SPECIES_STRONGEST_R2"])>=0.1 for x in q)
        for x in rows:
            if x["REGION_ID"]==reg["region_id"]:x["PRIMARY_SPECIES_SUPPORTING_SAME_LD_SIGNAL_N"]=n
    write_tsv(OUT/"PHASE14B_CORRECTED_SPECIES_LD_RESULTS.tsv",rows)
    write_tsv(OUT/"PHASE14B_HUNT_LOOKUP.tsv",hunt_rows)

    # Exact Phase11E-validated M1 projection decomposition at each disease lead.
    L=np.asarray([[float(r["Swedish_F1_loading_raw"]),float(r["Swedish_F2_loading_raw"])] for r in h2])
    R=np.loadtxt(BASE/"Swedish_R241_ldsc_raw_matrix.tsv",delimiter="\t",skiprows=1)
    D=np.sqrt(16017)*L; V=np.linalg.inv(D.T@np.linalg.solve(R,D)); A=V@D.T@np.linalg.inv(R); se_factor=np.sqrt(np.diag(V)); W=A/se_factor[:,None]
    multi=[]
    lead_indices=[]
    for reg in REGIONS:
        lead_rec=ref_rs[reg["lead"]]
        lead_indices.append(key_to_idx[canon(lead_rec[0],lead_rec[3],lead_rec[4],lead_rec[5])])
    lead_z_matrix=np.empty((len(manifest),len(REGIONS)),dtype=np.float64)
    for mi,m in enumerate(manifest):
        zz=np.load(m["swedish_matrix"],allow_pickle=False)["z"]
        lead_z_matrix[mi,:]=np.asarray(zz[lead_indices],dtype=np.float64)
    for reg in REGIONS:
        lead_rec=ref_rs[reg["lead"]]; lead_key=canon(lead_rec[0],lead_rec[3],lead_rec[4],lead_rec[5]); idx=key_to_idx[lead_key]
        zvec=lead_z_matrix[:,REGIONS.index(reg)]; contrib=W*zvec[None,:]; fz=contrib.sum(axis=1)
        for factor in range(2):
            absvals=np.abs(contrib[factor]); denom=absvals.sum(); order=np.argsort(-absvals)
            cancellation=1-abs(float(fz[factor]))/denom if denom else float("nan")
            for si,sp in enumerate(h2):
                if not (sp["Swedish_F1_top40"]=="TRUE" or sp["Swedish_F2_top40"]=="TRUE"):continue
                multi.append({"REGION_ID":reg["region_id"],"DISEASE":reg["disease"],"LEAD_SNP":reg["lead"],"FACTOR":f"F{factor+1}","FACTOR_LEAD_Z":float(fz[factor]),"SPECIES":sp["species"],"PRIMARY_TOP20_UNION":sp["Swedish_F1_top20"]=="TRUE" or sp["Swedish_F2_top20"]=="TRUE","SPECIES_LEAD_Z":float(zvec[si]),"FROZEN_LOADING":float(L[si,factor]),"EXACT_M1_SIGNED_CONTRIBUTION":float(contrib[factor,si]),"ABS_CONTRIBUTION_SHARE":float(absvals[si]/denom) if denom else float("nan"),"TOP1_ABS_SHARE":float(absvals[order[:1]].sum()/denom),"TOP5_ABS_SHARE":float(absvals[order[:5]].sum()/denom),"TOP20_ABS_SHARE":float(absvals[order[:20]].sum()/denom),"CANCELLATION_RATIO":cancellation,"OPERATOR_NOTE":"validated joint GLS operator with R_ldsc; not loading*Z"})
    write_tsv(OUT/"PHASE14B_MULTISPECIES_EFFECT_STRUCTURE.tsv",multi)

    # Full per-SNP coloc inputs for the six fixed candidate pairs.
    audit=[]
    for reg in REGIONS:
        sp=next(s for s in species if s["species"]==reg["candidate"])
        beta,se,z=candidate_arrays[sp["species"]]; dmap=disease_map[reg["region_id"]]; common=[]
        for key,d in dmap.items():
            idx=key_to_idx.get(key); rec=ref_key.get(key)
            if idx is None or rec is None:continue
            maf=maffun(rec)
            if not math.isfinite(maf) or maf<0.01:continue
            sk=str(keys[idx]);sa1=sk.split(":")[2];sa2=sk.split(":")[3]
            sb=float(beta[idx])
            if d["a1"]==sa1: aligned=sb
            elif d["a1"]==sa2:aligned=-sb
            else:continue
            common.append({"SNP_KEY":f"{key[0]}:{key[1]}:{key[2]}:{key[3]}","CHR":key[0],"BP":key[1],"DISEASE_A1":d["a1"],"DISEASE_A2":d["a2"],"DISEASE_BETA":d["beta"],"DISEASE_SE":d["se"],"DISEASE_P":d["p"],"DISEASE_INFO":d["info"],"DISEASE_NCASE":d["ncase"],"DISEASE_NCONTROL":d["ncontrol"],"SPECIES":sp["species"],"SPECIES_BETA_ALIGNED":aligned,"SPECIES_SE":float(se[idx]),"SPECIES_Z_ALIGNED":aligned/float(se[idx]),"SPECIES_P":pz(aligned/float(se[idx])),"REFERENCE_MAF_1000G_EUR":maf,"SPECIES_N_STUDY_LEVEL":16017,"SPECIES_INFO":"NA","SPECIES_SNP_SPECIFIC_N":"NA","MAF_PROVENANCE":"1000G EUR 489 reference; not original GWAS"})
        path=INPUTS/f"{reg['region_id']}__{sp['species'].replace(' ','_')}__coloc_input.tsv"
        write_tsv(path,common)
        audit.append({"REGION_ID":reg["region_id"],"DISEASE":reg["disease"],"CHR":reg["chrom"],"START":reg["start"],"END":reg["end"],"DISEASE_LEAD_SNP":reg["lead"],"PRESPECIFIED_SPECIES":reg["candidate"],"DISEASE_REGION_RAW_VARIANT_N":len(disease[reg["region_id"]]),"DENSE_SPECIES_REGION_VARIANT_N":len(region_indices[REGIONS.index(reg)]),"COLOC_COMMON_MAF_GE_0P01_SNP_N":len(common),"DISEASE_BETA_SE":"available","DISEASE_INFO":"available" if reg["disease"]=="SCZ" else "not supplied","DISEASE_N":"per-variant NCAS/NCON" if reg["disease"]=="SCZ" else "study-level cases/controls only","SPECIES_BETA_SE":"available","SPECIES_MAF":"1000G EUR reference sensitivity only","SPECIES_INFO":"not supplied","SPECIES_N":"study-level Swedish N=16017; SNP-specific N unavailable","GRCH_STATUS":"all analysis keys GRCh38","COLOC_INPUT":str(path),"FORMALITY":"coloc.abf sensitivity eligible; limitations must remain explicit" if len(common)>=100 else "insufficient common variants"})
    write_tsv(OUT/"PHASE14B_REGION_INPUT_AUDIT.tsv",audit)

    # Frozen six-cell peak context at candidate and high-LD proxy variants.
    beds=load_beds(); cell_rows=[]
    phase12={"SNCG_1":"SCZ Phase12A FDR-positive; Phase12B Swedish GLOBAL F2 Top5 negative conditional coefficient","SMC":"UC Phase12A FDR-positive","Macrophage":"CD Phase12A nominal","NaiveB":"CD Phase12A nominal","Early_Colonocyte":"UC Phase12A nominal","BNGA":"no frozen disease-coefficient highlight"}
    corrected={(x["REGION_ID"],x["SPECIES"]):x for x in rows}
    for reg in REGIONS:
        x=corrected[(reg["region_id"],reg["candidate"])]; ck=parse_key(x["SPECIES_STRONGEST_KEY"]); crec=ref_key.get(ck)
        proxies=[]
        for rec in ref_region_rows[reg["region_id"]]:
            rr=r2fun(crec,rec)
            if math.isfinite(rr) and rr>=0.8:proxies.append((rec,rr))
        for cell in beds:
            direct=in_peak(beds,cell,reg["chrom"],ck[1])
            hits=[(rec[2],rec[3],rr) for rec,rr in proxies if in_peak(beds,cell,reg["chrom"],rec[3])]
            cell_rows.append({"REGION_ID":reg["region_id"],"DISEASE":reg["disease"],"SPECIES":reg["candidate"],"CANDIDATE_KEY":x["SPECIES_STRONGEST_KEY"],"CELL":cell,"CANDIDATE_IN_PEAK":direct,"HIGH_LD_PROXY_R2_GE_0P8_N":len(proxies),"HIGH_LD_PROXY_IN_CELL_PEAK_N":len(hits),"HIGH_LD_PROXY_PEAK_EXAMPLES":";".join(f"{rid}:{bp}:r2={rr:.3f}" for rid,bp,rr in hits[:10]) or "NONE","CANDIDATE_GENE_CONTEXT":reg["gene"],"PHASE12_CONTEXT":phase12[cell],"EQTL_CAQTl_STATUS":"not available in frozen local resources; no molecular colocalization claimed","INTERPRETATION":"positional regulatory context only"})
    write_tsv(OUT/"PHASE14B_CELL_QTL_CONTEXT.tsv",cell_rows)
    for h in handles.values():h.close()
    (OUT/"checkpoints"/"phase14b_extraction.json").write_text(json.dumps({"status":"complete","fixed_regions":len(REGIONS),"species_rows":len(rows),"multispecies_rows":len(multi),"hunt_rows":len(hunt_rows),"cell_rows":len(cell_rows)},indent=2),encoding="utf-8")
    print(json.dumps({"status":"complete","corrected_rows":len(rows),"multispecies_rows":len(multi),"audit_rows":len(audit)}))


if __name__=="__main__":main()
