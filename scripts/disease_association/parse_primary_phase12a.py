#!/usr/bin/env python3
"""Parse completed Phase12A S-LDSC outputs, QC annotations, and render figures."""
import os
from pathlib import Path
import csv, gzip, json, math, re

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT=ROOT/'PHASE12A_CELL_DISEASE_SLD SC'
P9=ROOT/'FULL_SCALE_GFA_FINAL_EVALUATION'/'PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA'
REF=ROOT/'GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT'/'CELL_ANNOTATION_C1_EXECUTION'/'reference_download'/'GRCh38'/'baselineLD_v2.2_hm3'
CELLS=['brain__BNGA','brain__SNCG_1','brain__SMC','gut__Macrophage','gut__Early_Colonocyte','gut__NaiveB']
DISEASES=['SCZ','BD','AD','PD','CD','UC']
CELL_LABEL={'brain__BNGA':'BNGA','brain__SNCG_1':'SNCG_1','brain__SMC':'SMC','gut__Macrophage':'Macrophage','gut__Early_Colonocyte':'Early_Colonocyte','gut__NaiveB':'NaiveB'}
TISSUE={x: ('Brain' if x.startswith('brain') else 'Gut') for x in CELLS}

def read_tsv(p):
    with open(p,encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def f(x):
    try:return float(x)
    except:return float('nan')
def bh(rows, key='P'):
    good=sorted([(float(r[key]),i) for i,r in enumerate(rows) if r.get('STATUS')=='COMPLETED' and math.isfinite(float(r[key]))])
    n=len(good); q=[float('nan')]*len(rows); last=1.0
    for rank,(p,i) in reversed(list(enumerate(good,1))):
        last=min(last,p*n/rank);q[i]=last
    for i,x in enumerate(q): rows[i]['FDR_BH_30']='NA' if not math.isfinite(x) else '%.8g'%x
def main():
    manifest={r['annotation_id']:r for r in read_tsv(P9/'01_CELL_ATLAS_SOURCE_MANIFEST.tsv')}
    h2={}
    rows=[]
    # Standardized tau*: tau * SD(annotation) * M / h2; M is baseline HM3 M_5_50.
    m=sum(float(open(str(REF/('baselineLD.%d.l2.M_5_50'%c))).read().split()[0]) for c in range(1,23))
    for d in DISEASES:
      for cell in CELLS:
        q=dict(DISEASE=d,CELL_ID=cell,CELL_LABEL=CELL_LABEL[cell],TISSUE=TISSUE[cell],STATUS='COMPLETED',
               ANNOTATION_CLASS=('CELL_ACCESSIBLE' if cell.startswith('brain') else 'CELL_PREFERENTIAL_DAR'),
               S_LDSC_N_MODE=('NOT_RUN_NO_VERIFIED_N' if d=='PD' else 'SEE_PHASE12A_DISEASE_QC'))
        if d=='PD':
          q.update(STATUS='SKIPPED_NO_VERIFIED_N',REASON='Selected GP2 clinical-only PD input has no per-SNP or verified total N in frozen audit.')
          rows.append(q);continue
        res=OUT/'results'/(d+'__'+cell+'.results'); log=OUT/'logs'/(d+'__'+cell+'.h2.log')
        text=log.read_text(errors='replace'); mt=re.search(r'Total Observed scale h2: ([0-9.eE+-]+) \(([0-9.eE+-]+)\)',text)
        if not res.exists() or not mt:
          q.update(STATUS='FAILED_OR_MISSING',REASON='Missing .results or observed-scale h2 log line.');rows.append(q);continue
        h=float(mt.group(1)); se_h=float(mt.group(2)); h2[d]=(h,se_h)
        rr=read_tsv(res); hit=[x for x in rr if x['Category']=='L2_1']
        if len(hit)!=1: q.update(STATUS='FAILED_PARSE',REASON='Expected exactly one L2_1 category.');rows.append(q);continue
        x=hit[0]; p=math.erfc(abs(f(x['Coefficient_z-score']))/math.sqrt(2))
        prop=f(x['Prop._SNPs']); tau=f(x['Coefficient']); tau_se=f(x['Coefficient_std_error']); sd=math.sqrt(prop*(1-prop))
        q.update(H2_OBS=h,H2_OBS_SE=se_h,PROP_SNPS=prop,PROP_H2=f(x['Prop._h2']),PROP_H2_SE=f(x['Prop._h2_std_error']),ENRICHMENT=f(x['Enrichment']),ENRICHMENT_SE=f(x['Enrichment_std_error']),ENRICHMENT_P=f(x['Enrichment_p']),TAU=tau,TAU_SE=tau_se,TAU_STAR=tau*sd*m/h,TAU_STAR_SE=tau_se*sd*m/h,P=p,REASON='')
        rows.append(q)
    bh(rows)
    for q in rows:
      if q['STATUS']=='COMPLETED':
        q['TIER']='Tier 1' if float(q['FDR_BH_30'])<0.05 else ('Tier 2' if float(q['P'])<0.05 and ((q['DISEASE'] in ('SCZ','BD','AD') and q['TISSUE']=='Brain') or (q['DISEASE'] in ('CD','UC') and q['TISSUE']=='Gut')) else ('Tier 3' if float(q['P'])<0.10 else 'None'))
    fields=['DISEASE','CELL_ID','CELL_LABEL','TISSUE','ANNOTATION_CLASS','STATUS','S_LDSC_N_MODE','H2_OBS','H2_OBS_SE','PROP_SNPS','PROP_H2','PROP_H2_SE','ENRICHMENT','ENRICHMENT_SE','ENRICHMENT_P','TAU','TAU_SE','TAU_STAR','TAU_STAR_SE','P','FDR_BH_30','TIER','REASON']
    with open(OUT/'PHASE12A_36_PRIMARY_RESULTS.tsv','w',newline='',encoding='utf-8') as o:
      w=csv.DictWriter(o,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
    # annotation QC / overlap in exact HM3 baseline annotation space
    qc=[]; sets={}
    for cell in CELLS:
      s=set(); total=0
      for c in range(1,23):
        with gzip.open(OUT/'ldscores'/(cell+'.%d.annot.gz'%c),'rt') as h:
          next(h)
          for z in h:
            total+=1;v=z.rstrip().split('\t')[-1]
            if v=='1':s.add(z.rstrip().split('\t')[2])
      sets[cell]=s; mrow=manifest[cell]
      qc.append(dict(CELL_ID=cell,CELL_LABEL=CELL_LABEL[cell],TISSUE=TISSUE[cell],ANNOTATION_CLASS=('CELL_ACCESSIBLE' if cell.startswith('brain') else 'CELL_PREFERENTIAL_DAR'),PEAKS_ORIGINAL=mrow['original_assignment_peaks'],PEAKS_MERGED=mrow['merged_autosomal_intervals'],BP_MERGED=mrow['merged_autosomal_bp'],HM3_TOTAL=total,HM3_ANNOTATED=len(s),HM3_COVERAGE=len(s)/total,SOURCE_BED_SHA256=mrow['bed_sha256'],QC_STATUS='PASS'))
    with open(OUT/'PHASE12A_CELL_ANNOTATION_QC.tsv','w',newline='',encoding='utf-8') as o:
      w=csv.DictWriter(o,fieldnames=list(qc[0]),delimiter='\t');w.writeheader();w.writerows(qc)
    ov=[]
    for i,a in enumerate(CELLS):
      for b in CELLS[i+1:]:
        inter=len(sets[a]&sets[b]); union=len(sets[a]|sets[b]);ov.append(dict(CELL_A=a,CELL_B=b,HM3_SNP_INTERSECTION=inter,HM3_SNP_UNION=union,JACCARD=inter/union,OVERLAP_A=inter/len(sets[a]),OVERLAP_B=inter/len(sets[b])))
    with open(OUT/'PHASE12A_CELL_OVERLAP.tsv','w',newline='',encoding='utf-8') as o:
      w=csv.DictWriter(o,fieldnames=list(ov[0]),delimiter='\t');w.writeheader();w.writerows(ov)
    # exact figures: quantitative grid, primary tau* as hero evidence; no result selection.
    # Imported only after tabular outputs so a plotting-runtime fault never
    # blocks statistical parsing.
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'svg.fonttype':'none','pdf.fonttype':42,'font.size':7,'axes.spines.right':False,'axes.spines.top':False})
    done=[x for x in rows if x['STATUS']=='COMPLETED']; mat=np.full((6,6),np.nan)
    for x in done:mat[DISEASES.index(x['DISEASE']),CELLS.index(x['CELL_ID'])]=float(x['TAU_STAR'])
    fig,ax=plt.subplots(figsize=(7.2,3.8));im=ax.imshow(mat,cmap='RdBu_r',aspect='auto');fig.colorbar(im,ax=ax,label='Standardized effect (tau*)');ax.set_xticks(range(6),[CELL_LABEL[x] for x in CELLS],rotation=30,ha='right');ax.set_yticks(range(6),DISEASES);ax.set_title('Primary S-LDSC: baselineLD-adjusted cell-annotation tau*');
    for i in range(6):
      for j in range(6):
       if np.isfinite(mat[i,j]):ax.text(j,i,('*' if float([x for x in done if x['DISEASE']==DISEASES[i] and x['CELL_ID']==CELLS[j]][0]['FDR_BH_30'])<.05 else ''),ha='center',va='center',fontsize=9)
    for ext in ('svg','pdf','png'):fig.savefig(OUT/('PHASE12A_FIGURE_1_CELL_DISEASE_HEATMAP.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(7.2,3.5),sharey=True)
    for ax,tissue in zip(axs,['Brain','Gut']):
      sub=[x for x in done if x['TISSUE']==tissue];labs=[x['DISEASE']+'\n'+x['CELL_LABEL'] for x in sub]; vals=[float(x['ENRICHMENT']) for x in sub];errs=[float(x['ENRICHMENT_SE']) for x in sub];colors=['#4878A8' if tissue=='Brain' else '#59A14F']*len(sub);ax.bar(range(len(sub)),vals,yerr=errs,color=colors,ecolor='#333333',capsize=2);ax.axhline(1,color='#666',lw=.8);ax.set_title(tissue+' annotations');ax.set_xticks(range(len(sub)),labs,rotation=90,fontsize=5);ax.set_ylabel('Heritability enrichment')
    fig.suptitle('Marginal enrichment by tissue context (baselineLD-adjusted coefficient is primary)');fig.tight_layout()
    for ext in ('svg','pdf','png'):fig.savefig(OUT/('PHASE12A_FIGURE_2_TISSUE_SPECIFIC_ENRICHMENT.'+ext),dpi=600,bbox_inches='tight')
if __name__=='__main__':main()
