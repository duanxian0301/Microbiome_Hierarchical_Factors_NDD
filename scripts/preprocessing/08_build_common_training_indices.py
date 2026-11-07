"""Freeze the actual cross-cohort, LD-controlled training order and indices.

This is identity bookkeeping only: Swedish establishes the analysis order;
HUNT is explicitly re-indexed to that order.  The resulting files are used by
all downstream matrix assembly, preventing an accidental positional join.
"""
import csv
import hashlib
import os
from pathlib import Path

import numpy as np

ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
OUT = ROOT / "FULL_SCALE_GFA_FINAL_EVALUATION"
manifest = list(csv.DictReader(open(OUT / "inputs_241_manifest.tsv", encoding="utf-8"), delimiter="\t"))
with np.load(manifest[0]["swedish_matrix"]) as x:
    swedish = x["snp"]
with np.load(manifest[0]["hunt_matrix"]) as x:
    hunt = x["snp"]

training_rows = list(csv.DictReader(open(OUT / "genomewide_ld_controlled_training_manifest.tsv", encoding="utf-8"), delimiter="\t"))
rows_by_key = {row["canonical_key"]: row for row in training_rows}
if len(rows_by_key) != len(training_rows):
    raise RuntimeError("Duplicate canonical key in LD-controlled manifest")

# Dense GWAS keys preserve REF/ALT order; the resource key sorts alleles.  Thus
# a direct string join is deliberately rejected and identity is reconstructed
# by chromosome, position, and unordered allele pair.
def resource_key(dense_key):
    chrom, pos, a1, a2 = dense_key.split(":")
    a1, a2 = sorted((a1, a2))
    return f"{chrom}:{pos}:{a1}:{a2}"

swedish_by_resource = {}
for i, key in enumerate(swedish):
    ck = resource_key(key)
    if ck in rows_by_key:
        if ck in swedish_by_resource:
            raise RuntimeError(f"Duplicate dense Swedish representation for {ck}")
        swedish_by_resource[ck] = i
if len(swedish_by_resource) != len(training_rows):
    raise RuntimeError(f"Swedish LD resource mismatch: selected={len(swedish_by_resource)}, manifest={len(training_rows)}")

# Preserve dense Swedish order; re-index HUNT by the exact dense key, which
# had already been independently verified to have identical membership.
swedish_index = np.array(sorted(swedish_by_resource.values()), dtype=np.int32)
keys = swedish[swedish_index]
needed_dense_keys = set(keys)
hunt_by_dense = {}
for i, key in enumerate(hunt):
    if key in needed_dense_keys:
        if key in hunt_by_dense:
            raise RuntimeError(f"Duplicate dense HUNT key {key}")
        hunt_by_dense[key] = i
hunt_index = np.fromiter((hunt_by_dense.get(key, -1) for key in keys), dtype=np.int32)
if (hunt_index < 0).any():
    raise RuntimeError(f"HUNT dense-key mismatch: missing={(hunt_index < 0).sum()}")

# Manifest follows fixed Swedish dense order, not chromosome file traversal.
ordered_rows = [rows_by_key[resource_key(key)] for key in keys]

np.save(OUT / "Swedish_common_training_index.npy", swedish_index)
np.save(OUT / "HUNT_common_training_index.npy", hunt_index)
np.save(OUT / "common_training_dense_key.npy", keys)
with open(OUT / "common_genomewide_ld_controlled_training_manifest.tsv", "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(training_rows[0]), delimiter="\t")
    writer.writeheader(); writer.writerows(ordered_rows)

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

report = OUT / "common_training_identity_report.txt"
report.write_text(
    "swedish_dense_variants=%d\n"
    "hunt_dense_variants=%d\n"
    "training_variants=%d\n"
    "swedish_order_is_hunt_order=%s\n"
    "all_hunt_indices_resolved=%s\n"
    "dense_key_sha256=%s\n"
    "swedish_index_sha256=%s\n"
    "hunt_index_sha256=%s\n"
    % (len(swedish), len(hunt), len(keys), bool(np.array_equal(swedish, hunt)), bool((hunt_index >= 0).all()),
       sha(OUT / "common_training_dense_key.npy"), sha(OUT / "Swedish_common_training_index.npy"), sha(OUT / "HUNT_common_training_index.npy")),
    encoding="utf-8",
)
print(report.read_text(encoding="utf-8"), end="")
