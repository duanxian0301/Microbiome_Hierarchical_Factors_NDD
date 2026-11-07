#!/usr/bin/env python3
"""Reproducible Phase 12A primary S-LDSC execution.

Uses the frozen Phase 9 BED annotations without changing intervals, an existing
GRCh38 1000G-EUR/baselineLD-v2.2/HM3 bundle, and Phase11C harmonized disease
statistics.  PD is intentionally not munged: the selected clinical-only source
does not supply a defensible sample size in the frozen input/audit.
"""
from __future__ import annotations
import csv, gzip, json, math, os, re, shutil, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import os
from pathlib import Path

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT = ROOT / "PHASE12A_CELL_DISEASE_SLD SC"
P9 = ROOT / "FULL_SCALE_GFA_FINAL_EVALUATION" / "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA"
HARM = ROOT / "PHASE11_DISEASE_BIOLOGICAL_INTEGRATION" / "02_DISEASE_RESOURCES" / "PHASE11C_HARMONIZED"
LDSCROOT = ROOT / "GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT" / "CELL_ANNOTATION_C1_EXECUTION"
REF = LDSCROOT / "reference_download"
PY27 = Path(os.environ["MHFNDD_LDSC_PYTHON"]).resolve()
LDSC = LDSCROOT / "software" / "ldsc" / "ldsc.py"
MUNGE = LDSCROOT / "software" / "ldsc" / "munge_sumstats.py"
BUILDER = LDSCROOT / "scripts" / "build_binary_peak_annotations.py"
BFILE = REF / "GRCh38" / "plink_files"
BASELINE = REF / "GRCh38" / "baselineLD_v2.2_hm3" / "baselineLD."
WEIGHTS = REF / "GRCh38" / "weights" / "weights.hm3_noMHC."
HM3 = REF / "weights_hm3_noMHC" / "w_hm3.snplist"

PRIMARY = ["brain__BNGA", "brain__SNCG_1", "brain__SMC", "gut__Macrophage", "gut__Early_Colonocyte", "gut__NaiveB"]
# Per disease fixed effective N when per-SNP sample size is absent.  Values are
# determined solely from the audited Phase11B case/control counts.
FIXED_CASE_CONTROL = {"CD": (12194, 28072), "UC": (12366, 33609)}
ALL_DISEASES = ["SCZ", "BD", "AD", "PD", "CD", "UC"]

def run(command, log):
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "ab") as h:
        h.write(("\nCOMMAND: " + " ".join(map(str, command)) + "\n").encode())
        p = subprocess.run([str(x) for x in command], stdout=h, stderr=subprocess.STDOUT, env=os.environ.copy())
    if p.returncode:
        raise RuntimeError("failed (%s): %s" % (p.returncode, command))

def read_tsv(path):
    with open(path, newline='', encoding='utf-8-sig') as h:
        return list(csv.DictReader(h, delimiter='\t'))

def primary_metadata():
    rows = {r["annotation_id"]: r for r in read_tsv(P9 / "01_CELL_ATLAS_SOURCE_MANIFEST.tsv")}
    return [rows[x] for x in PRIMARY]

def build_annots():
    ann = OUT / "annotations"
    for row in primary_metadata():
        name, bed = row["annotation_id"], Path(row["bed_path"])
        if not bed.exists(): raise FileNotFoundError(bed)
        if not (ann / (name + ".1.annot.gz")).exists():
            run([sys.executable, BUILDER, "--bed", bed, "--bim-dir", BFILE, "--out-dir", ann, "--prefix", name], OUT / "logs" / (name + ".annotation.log"))

def make_ldscores():
    dest = OUT / "ldscores"
    dest.mkdir(parents=True, exist_ok=True)
    jobs=[]
    for name in PRIMARY:
        for chrom in range(1, 23):
            pre = dest / (name + "." + str(chrom))
            if (Path(str(pre) + ".l2.ldscore.gz")).exists(): continue
            jobs.append((name, chrom, pre))
    def one(job):
        name, chrom, pre = job
        run([PY27, LDSC, "--l2", "--bfile", BFILE / ("1000G.EUR.hg38.%d" % chrom), "--ld-wind-cm", "1", "--annot", OUT / "annotations" / (name + ".%d.annot.gz" % chrom), "--thin-annot", "--print-snps", REF / "regression_snps_grch38" / ("hm3_noMHC.%d.snplist" % chrom), "--out", pre], OUT / "logs" / (name + ".%d.l2.log" % chrom))
    # Four bounded independent chromosome jobs make use of the local machine
    # without oversubscribing it; each completed output is a durable checkpoint.
    with ThreadPoolExecutor(max_workers=4) as ex:
        for f in as_completed([ex.submit(one, x) for x in jobs]):
            f.result()

def bim_map():
    # Mapping disease CHR:BP:A1:A2 keys to LDSC rsids. Require a unique unordered allele match.
    d = defaultdict(list)
    for chrom in range(1, 23):
        with open(BFILE / ("1000G.EUR.hg38.%d.bim" % chrom)) as h:
            for z in h:
                c, rs, cm, bp, a1, a2 = z.split()[:6]
                d[(c, bp, tuple(sorted((a1.upper(), a2.upper()))))].append(rs)
    return {k:v[0] for k,v in d.items() if len(v)==1}

def neff(ncase, ncontrol): return 4.0 / (1.0/ncase + 1.0/ncontrol)

def prepare_sumstats():
    outdir = OUT / "sumstats_input"; outdir.mkdir(parents=True, exist_ok=True)
    # On a resumed run, retain the already audited/munged five-disease set and
    # avoid rebuilding the multi-million-SNP BIM lookup solely to rediscover it.
    ready = [d for d in ALL_DISEASES if (outdir/(d+'.tsv')).exists() and (OUT/'munged'/(d+'.sumstats.gz')).exists()]
    if ready:
        return ready
    bmap = bim_map()
    qc=[]
    for disease in ALL_DISEASES:
        source = HARM / (disease + "_HM3_GRCh38.tsv.gz")
        if disease == "PD":
            qc.append(dict(disease=disease, status="SKIPPED_NO_VERIFIED_N", input_rows="571558", mapped_rows="NA", n_mode="NO_INFERENTIAL_N", N="NA", reason="Selected GP2 clinical-only input has no per-SNP or verified total case/control N in frozen audit.")); continue
        total=mapped=invalid_n=0; ns=[]; target=outdir/(disease+".tsv")
        with gzip.open(source, 'rt') as src, open(target, 'w', newline='') as dst:
            r=csv.DictReader(src, delimiter='\t'); w=csv.DictWriter(dst, fieldnames=['SNP','A1','A2','BETA','P','N'], delimiter='\t'); w.writeheader()
            for x in r:
                total+=1; key=(x['CHR'],x['BP'],tuple(sorted((x['A1'].upper(),x['A2'].upper())))); rs=bmap.get(key)
                if not rs: continue
                try:
                    if disease in FIXED_CASE_CONTROL: n=neff(*FIXED_CASE_CONTROL[disease]); mode='fixed_case_control_Neff'
                    else:
                        ncase=float(x['N_CASE']); nctrl=float(x['N_CONTROL']); n=neff(ncase,nctrl); mode='per_SNP_Neff'
                    if not (n>0 and math.isfinite(n)): raise ValueError()
                    beta=float(x['BETA']); p=float(x['P'])
                    if not (math.isfinite(beta) and 0 < p <= 1): raise ValueError()
                except (ValueError, KeyError): invalid_n+=1; continue
                w.writerow(dict(SNP=rs,A1=x['A1'],A2=x['A2'],BETA=beta,P=p,N='%.6f'%n)); mapped+=1; ns.append(n)
        qc.append(dict(disease=disease,status='READY',input_rows=total,mapped_rows=mapped,n_mode=mode,N_median=sorted(ns)[len(ns)//2] if ns else 'NA',invalid_or_unmapped=total-mapped,reason=''))
    with open(OUT/'PHASE12A_DISEASE_QC.tsv','w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=sorted({k for q in qc for k in q}),delimiter='\t',extrasaction='ignore'); w.writeheader(); w.writerows(qc)
    return [x['disease'] for x in qc if x['status']=='READY']

def munge(diseases):
    dest=OUT/'munged'; dest.mkdir(parents=True,exist_ok=True)
    for d in diseases:
        pre=dest/d
        if Path(str(pre)+'.sumstats.gz').exists(): continue
        run([PY27,MUNGE,'--sumstats',OUT/'sumstats_input'/(d+'.tsv'),'--snp','SNP','--a1','A1','--a2','A2','--signed-sumstats','BETA,0','--p','P','--N-col','N','--merge-alleles',HM3,'--out',pre],OUT/'logs'/(d+'.munge.log'))

def h2(diseases):
    dest=OUT/'results'; dest.mkdir(parents=True,exist_ok=True)
    # LDSC --overlap-annot reconstructs the overlap matrix from annotation
    # files resolved next to each custom LD-score prefix.  It requires the BIM
    # key columns as well as the binary column (the thin annot used for --l2
    # deliberately has only ANNOT).  Recreate keyed immutable copies here.
    for cell in PRIMARY:
        for chrom in range(1, 23):
            src=OUT/'annotations'/(cell+'.%d.annot.gz' % chrom)
            dst=OUT/'ldscores'/(cell+'.%d.annot.gz' % chrom)
            # baselineLD_v2.2_hm3 annotations contain only its HM3 regression
            # rows, whereas the bfile contains all 1000G variants.  Match the
            # custom values to that exact baseline annotation order; otherwise
            # LDSC cannot construct the overlap matrix.
            base=REF/'GRCh38'/'baselineLD_v2.2_hm3'/("baselineLD.%d.annot.gz" % chrom)
            base_rows=[]; wanted=set()
            with gzip.open(base, 'rt') as inp:
                header=next(inp).rstrip('\n').split('\t')
                for line in inp:
                    f=line.rstrip('\n').split('\t'); base_rows.append(f); wanted.add(f[2])
            values={}
            with gzip.open(src, 'rt') as inp, open(BFILE/("1000G.EUR.hg38.%d.bim" % chrom)) as bim:
                next(inp)
                for ann, line in zip(inp, bim):
                    f=line.split()
                    if f[1] in wanted: values[f[1]]=ann.strip()
            if len(values) != len(wanted): raise RuntimeError('custom/base HM3 mismatch %s chr%s: %s/%s' % (cell, chrom, len(values), len(wanted)))
            with gzip.open(dst, 'wt') as out:
                out.write('CHR\tBP\tSNP\tCM\tANNOT\n')
                for f in base_rows: out.write('%s\t%s\t%s\t%s\t%s\n' % (f[0], f[1], f[2], f[3], values[f[2]]))
    for d in diseases:
        for cell in PRIMARY:
            pre=dest/(d+'__'+cell)
            if Path(str(pre)+'.results').exists(): continue
            run([PY27,LDSC,'--h2',OUT/'munged'/(d+'.sumstats.gz'),'--ref-ld-chr',str(BASELINE)+','+str(OUT/'ldscores'/cell)+'.','--w-ld-chr',WEIGHTS,'--frqfile-chr',str(BFILE/'1000G.EUR.hg38.'),'--overlap-annot','--print-coefficients','--out',pre],OUT/'logs'/(d+'__'+cell+'.h2.log'))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    build_annots(); make_ldscores(); ds=prepare_sumstats(); munge(ds); h2(ds)
    (OUT/'RUN_PRIMARY_COMPLETE.flag').write_text('Complete; parse with parse_phase12a_results.py\n')

if __name__=='__main__': main()
