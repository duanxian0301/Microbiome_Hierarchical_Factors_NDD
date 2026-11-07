"""Build exact-order 241-trait inputs for the original official GFA R_ldsc code."""
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT = ROOT / "FULL_SCALE_GFA_FINAL_EVALUATION"
INP = OUT / "r241_inputs"
INP.mkdir(exist_ok=True)
rows = list(csv.DictReader(open(OUT / "inputs_241_manifest.tsv", encoding="utf-8"), delimiter="\t"))
ref = np.load(ROOT / "GFA_NATIVE_R_REESTIMATION_AND_VALIDATION" / "inputs" / "official_hm3_nomHC_ldscore_canonical_keys.npz")
refkey = ref["canonical_key"]
refl2 = ref["L2"].astype(np.float64)

with np.load(rows[0]["swedish_matrix"]) as x:
    swedish_all = x["snp"]
with np.load(rows[0]["hunt_matrix"]) as x:
    hunt_all = x["snp"]

# The LDSC reference was historically constructed in the dense key orientation.
# Exact key matching is intentional here; the original GFA R_ldsc routine needs
# L2 for this exact allele-verified reference identity, not a positional join.
swedish_index = np.flatnonzero(np.isin(swedish_all, refkey, assume_unique=False)).astype(np.int32)
keys = swedish_all[swedish_index]
needed = set(keys)
hunt_lookup = {}
for i, key in enumerate(hunt_all):
    if key in needed:
        if key in hunt_lookup:
            raise RuntimeError(f"Duplicate HUNT dense key: {key}")
        hunt_lookup[key] = i
hunt_index = np.fromiter((hunt_lookup.get(k, -1) for k in keys), dtype=np.int32)
if (hunt_index < 0).any():
    raise RuntimeError(f"HUNT LDSC input missing {(hunt_index < 0).sum()} Swedish-reference variants")

order = np.argsort(refkey)
where = np.searchsorted(refkey[order], keys)
if not np.array_equal(refkey[order][where], keys):
    raise RuntimeError("LDSC L2 reference alignment failure")
l2 = refl2[order][where]
if not (np.isfinite(l2).all() and (l2 > 0).all()):
    raise RuntimeError("Invalid LDSC L2")
np.save(INP / "Swedish_ldsc_index.npy", swedish_index)
np.save(INP / "HUNT_ldsc_index.npy", hunt_index)
np.save(INP / "ldsc_dense_key.npy", keys)
np.ascontiguousarray(l2, dtype=np.float64).tofile(INP / "hm3_nomHC_L2.float64.bin")

audits = []
for cohort, matrix_field, index in (("Swedish", "swedish_matrix", swedish_index), ("HUNT", "hunt_matrix", hunt_index)):
    # R reads column-major data; writing transpose C-contiguous produces exactly that layout.
    z = np.empty((len(keys), len(rows)), dtype=np.float64)
    for j, row in enumerate(rows):
        with np.load(row[matrix_field]) as x:
            zz = x["z"][index].astype(np.float64)
            bb = x["beta"][index].astype(np.float64)
            se = x["se"][index].astype(np.float64)
        valid = np.isfinite(zz) & np.isfinite(bb) & np.isfinite(se) & (se > 0)
        if not valid.all():
            raise RuntimeError(f"{cohort} trait {j+1}: invalid LDSC input")
        delta = np.abs(zz - bb / se)
        audits.append({"cohort": cohort, "trait_order": row["trait_order"], "accession": row["swedish_accession" if cohort == "Swedish" else "hunt_accession"],
                       "ldsc_snps": len(keys), "finite": int(valid.sum()), "z_beta_se_max_abs_diff": format(float(delta.max()), ".8g")})
        z[:, j] = zz
    np.ascontiguousarray(z.T, dtype=np.float64).tofile(INP / f"{cohort}_Z_ldsc.float64.bin")
    del z

with open(INP / "R241_INPUT_TRAIT_AUDIT.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(audits[0]), delimiter="\t"); w.writeheader(); w.writerows(audits)
MROOT = ROOT / "GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT" / "CELL_ANNOTATION_C1_EXECUTION" / "reference_download" / "GRCh38" / "weights"
M = sum(int((MROOT / f"weights.hm3_noMHC.{c}.l2.M").read_text().strip()) for c in range(1, 23))
meta = {"J_common_allele_verified": int(len(keys)), "K_traits": len(rows), "ld_size_sum_M_official_weight_reference": M,
        "trait_order": [r["trait_order"] for r in rows], "layout": "float64 trait-contiguous for R matrix(readBin,nrow=J,ncol=K)"}
(INP / "R241_input_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()
with open(INP / "R241_INPUT_SHA256.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(("file", "sha256"))
    for p in (INP / "Swedish_Z_ldsc.float64.bin", INP / "HUNT_Z_ldsc.float64.bin", INP / "hm3_nomHC_L2.float64.bin", INP / "ldsc_dense_key.npy"):
        w.writerow((p.name, sha(p)))
print(f"ldsc_snps={len(keys)} traits={len(rows)} Swedish/HUNT_order_reconciled=true")
