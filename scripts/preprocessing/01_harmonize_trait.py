"""Checkpointed single-trait dense harmonization; safe to run independently."""
import argparse, gzip, hashlib, json, math
import os
from pathlib import Path
import numpy as np
import pandas as pd
from pyliftover import LiftOver

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT=ROOT/'CELL_ANNOTATION_FACTOR_LOCALIZATION'/'DENSE_HARMONIZED_GWAS_EXECUTION'
REF=ROOT/'GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT'/'CELL_ANNOTATION_C1_EXECUTION'/'reference_download'/'GRCh38'/'plink_files'
MAN = Path(os.environ["MHFNDD_TRAIT_MANIFEST"]).resolve()
IDS = Path(os.environ["MHFNDD_COMMON_SNP_IDS"]).resolve()
CHAIN = Path(os.environ["MHFNDD_HG19_TO_HG38_CHAIN"]).resolve()
def reference():
    targets=set(IDS.read_text(encoding='utf-8').splitlines()); byrs={}; bycoord={}
    for chrom in range(1,23):
        with open(REF/f'1000G.EUR.hg38.{chrom}.bim',encoding='utf-8') as h:
            for line in h:
                z=line.split()
                if len(z)!=6 or not z[0].isdigit() or z[1] not in targets: continue
                rec=(int(z[0]),int(float(z[3])),z[4].upper(),z[5].upper())
                byrs.setdefault(z[1],rec); bycoord.setdefault((rec[0],rec[1],frozenset((rec[2],rec[3]))),z[1])
    return byrs,bycoord
def pal(a,b): return {a,b} in ({'A','T'},{'C','G'})
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--cohort',choices=['Swedish','HUNT'],required=True); ap.add_argument('--index',type=int,required=True); a=ap.parse_args()
    m=pd.read_csv(MAN,sep='\t').drop_duplicates(['swedish','hunt']).reset_index(drop=True); row=m.iloc[a.index]
    accession=row.swedish if a.cohort=='Swedish' else row.hunt; path=Path(row.swedish_path if a.cohort=='Swedish' else row.hunt_path)
    od=OUT/'blocks'/a.cohort; od.mkdir(parents=True,exist_ok=True); npz=od/f'{accession}.npz'; qc=od/f'{accession}.qc.json'
    if npz.exists() and qc.exists(): return
    byrs,bycoord=reference(); lo=LiftOver(str(CHAIN)) if a.cohort=='Swedish' else None
    snpi=[]; bet=[]; seout=[]; zval=[]; counts={k:0 for k in ['records','reference_candidate','coord_ok','allele_match','palindromic_removed','invalid_statistics','retained','flipped','unmapped_liftover']}
    with gzip.open(path,'rt',encoding='utf-8',newline='') as h:
        ix={n:i for i,n in enumerate(h.readline().rstrip('\r\n').split('\t'))}; idfield='rs_id' if 'rs_id' in ix else None
        needed={'chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error'}
        if not needed.issubset(ix): raise ValueError(f'{path}: missing {needed-set(ix)}')
        for line in h:
            counts['records']+=1; q=line.rstrip('\r\n').split('\t')
            try:
                c=int(q[ix['chromosome']]); bp=int(float(q[ix['base_pair_location']])); ea=q[ix['effect_allele']].upper(); oa=q[ix['other_allele']].upper(); b=float(q[ix['beta']]); s=float(q[ix['standard_error']])
            except: continue
            if idfield:
                rs=q[ix[idfield]]; rec=byrs.get(rs)
            else:
                rs=bycoord.get((c,bp,frozenset((ea,oa)))); rec=byrs.get(rs) if rs else None
            if rec is None: continue
            counts['reference_candidate']+=1
            if a.cohort=='Swedish':
                hit=lo.convert_coordinate('chr'+str(c),bp-1)
                ok=bool(hit) and hit[0][0]=='chr'+str(rec[0]) and hit[0][1]+1==rec[1]
                if not hit: counts['unmapped_liftover']+=1
            else: ok=(c==rec[0] and bp==rec[1])
            if not ok: continue
            counts['coord_ok']+=1
            if {ea,oa}!={rec[2],rec[3]}: continue
            counts['allele_match']+=1
            if pal(ea,oa): counts['palindromic_removed']+=1; continue
            if not (math.isfinite(b) and math.isfinite(s) and s>0): counts['invalid_statistics']+=1; continue
            if ea!=rec[2]: b=-b; counts['flipped']+=1
            snpi.append(rs); bet.append(b); seout.append(s); zval.append(b/s); counts['retained']+=1
    np.savez_compressed(npz,snp=np.asarray(snpi,dtype='U32'),beta=np.asarray(bet,dtype='f4'),se=np.asarray(seout,dtype='f4'),z=np.asarray(zval,dtype='f4'))
    counts.update(cohort=a.cohort,accession=accession,species=row.species,path=str(path),reference='GRCh38_1000G_EUR',coordinate_mode='hg19ToHg38' if a.cohort=='Swedish' else 'direct_GRCh38',chain=str(CHAIN) if a.cohort=='Swedish' else None)
    qc.write_text(json.dumps(counts,indent=2),encoding='utf-8')
if __name__=='__main__': main()
