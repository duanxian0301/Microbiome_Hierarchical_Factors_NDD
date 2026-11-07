#!/usr/bin/env python3
"""Phase17 exact Phase12A-framework extension for ET and ITL6_2_2."""
from __future__ import annotations
import csv, gzip, math, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT=ROOT/"PHASE17_EXPLORATORY_CELL_DISEASE_EXTENSION"
P9=ROOT/"FULL_SCALE_GFA_FINAL_EVALUATION"/"PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"
P12=ROOT/"PHASE12A_CELL_DISEASE_SLD SC"
LDSCROOT=ROOT/"GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT"/"CELL_ANNOTATION_C1_EXECUTION"
REF=LDSCROOT/"reference_download"
PY27 = Path(os.environ["MHFNDD_LDSC_PYTHON"]).resolve()
LDSC=LDSCROOT/"software"/"ldsc"/"ldsc.py"
BUILDER=LDSCROOT/"scripts"/"build_binary_peak_annotations.py"
BFILE=REF/"GRCh38"/"plink_files"
BASELINE=REF/"GRCh38"/"baselineLD_v2.2_hm3"/"baselineLD."
WEIGHTS=REF/"GRCh38"/"weights"/"weights.hm3_noMHC."
REG=REF/"regression_snps_grch38"
CELLS=["brain__ET","brain__ITL6_2_2"]
DISEASES=["SCZ","BD","AD","CD","UC"]

def tsv(path):
    with open(path,encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h,delimiter="\t"))
def run(cmd,log):
    log.parent.mkdir(parents=True,exist_ok=True)
    with open(log,"ab") as h:
        h.write(("\nCOMMAND: "+" ".join(map(str,cmd))+"\n").encode())
        p=subprocess.run([str(x) for x in cmd],stdout=h,stderr=subprocess.STDOUT,
            env=os.environ.copy())
    if p.returncode: raise RuntimeError("failed %s; see %s"%(cmd,log))

def manifest(): return {x["annotation_id"]:x for x in tsv(P9/"01_CELL_ATLAS_SOURCE_MANIFEST.tsv")}

def build():
    meta=manifest()
    for cell in CELLS:
        row=meta[cell]; bed=Path(row["bed_path"])
        if not bed.exists(): raise FileNotFoundError(bed)
        if not (OUT/"annotations"/(cell+".1.annot.gz")).exists():
            run([sys.executable,BUILDER,"--bed",bed,"--bim-dir",BFILE,"--out-dir",OUT/"annotations","--prefix",cell],OUT/"logs"/(cell+".annotation.log"))

def l2():
    jobs=[]
    for cell in CELLS:
      for chrom in range(1,23):
        pre=OUT/"ldscores"/(cell+"."+str(chrom))
        if not Path(str(pre)+".l2.ldscore.gz").exists():jobs.append((cell,chrom,pre))
    def one(j):
      cell,chrom,pre=j
      run([PY27,LDSC,"--l2","--bfile",BFILE/("1000G.EUR.hg38.%d"%chrom),"--ld-wind-cm","1","--annot",OUT/"annotations"/(cell+".%d.annot.gz"%chrom),"--thin-annot","--print-snps",REG/("hm3_noMHC.%d.snplist"%chrom),"--out",pre],OUT/"logs"/(cell+".%d.l2.log"%chrom))
    with ThreadPoolExecutor(max_workers=4) as ex:
      for f in as_completed([ex.submit(one,j) for j in jobs]):f.result()

def keyed_annots():
    for cell in CELLS:
      for chrom in range(1,23):
        src=OUT/"annotations"/(cell+".%d.annot.gz"%chrom)
        dst=OUT/"ldscores"/(cell+".%d.annot.gz"%chrom)
        base=REF/"GRCh38"/"baselineLD_v2.2_hm3"/("baselineLD.%d.annot.gz"%chrom)
        rows=[];wanted=set()
        with gzip.open(base,"rt") as h:
          header=next(h)
          for line in h:
            f=line.rstrip("\n").split("\t");rows.append(f);wanted.add(f[2])
        vals={}
        with gzip.open(src,"rt") as a,open(BFILE/("1000G.EUR.hg38.%d.bim"%chrom)) as b:
          next(a)
          for av,bv in zip(a,b):
            bf=bv.split()
            if bf[1] in wanted:vals[bf[1]]=av.strip()
        if len(vals)!=len(wanted):raise RuntimeError("baseline/custom mismatch %s chr%d"%(cell,chrom))
        with gzip.open(dst,"wt") as h:
          h.write("CHR\tBP\tSNP\tCM\tANNOT\n")
          for f in rows:h.write("%s\t%s\t%s\t%s\t%s\n"%(f[0],f[1],f[2],f[3],vals[f[2]]))

def regress():
    for d in DISEASES:
      ss=P12/"munged"/(d+".sumstats.gz")
      if not ss.exists():raise FileNotFoundError(ss)
      for cell in CELLS:
        pre=OUT/"results"/(d+"__"+cell)
        if Path(str(pre)+".results").exists():continue
        run([PY27,LDSC,"--h2",ss,"--ref-ld-chr",str(BASELINE)+","+str(OUT/"ldscores"/cell)+".","--w-ld-chr",WEIGHTS,"--frqfile-chr",str(BFILE/"1000G.EUR.hg38."),"--overlap-annot","--print-coefficients","--out",pre],OUT/"logs"/(d+"__"+cell+".h2.log"))

def main():
    build();l2();keyed_annots();regress();(OUT/"RUN_COMPLETE.flag").write_text("Phase17 2x5 S-LDSC complete\n",encoding="utf-8")
if __name__=="__main__":main()
