"""Run fixed two-step PLINK2 sort then r2 pruning for one annotation/chromosome."""
import argparse, subprocess
import os
from pathlib import Path

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT=ROOT/'CELL_ANNOTATION_FACTOR_LOCALIZATION'/'DENSE_HARMONIZED_GWAS_EXECUTION'
REF=ROOT/'GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT'/'CELL_ANNOTATION_C1_EXECUTION'/'reference_download'/'GRCh38'/'plink_files'
P=OUT/'tools'/'plink2_20261007'/'plink2.exe'
ap=argparse.ArgumentParser(); ap.add_argument('--annotation',required=True); ap.add_argument('--chromosome',type=int,required=True); a=ap.parse_args()
wd=OUT/'ld_pruning'/'runs'/a.annotation; wd.mkdir(parents=True,exist_ok=True)
pre=wd/f'chr{a.chromosome}_sorted'; final=wd/f'chr{a.chromosome}_r2_02'
if not final.with_suffix('.prune.in').exists():
    if not pre.with_suffix('.pgen').exists():
        subprocess.run([str(P),'--bfile',str(REF/f'1000G.EUR.hg38.{a.chromosome}'),'--extract',str(OUT/'ld_pruning'/'extract_lists'/a.annotation/f'chr{a.chromosome}.snps'),'--make-pgen','--sort-vars','--out',str(pre),'--threads','4'],check=True)
    subprocess.run([str(P),'--pfile',str(pre),'--indep-pairwise','250kb','1','0.2','--out',str(final),'--threads','4'],check=True)
