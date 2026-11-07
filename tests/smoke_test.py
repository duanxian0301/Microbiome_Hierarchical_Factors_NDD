#!/usr/bin/env python3
"""Lightweight public-package integrity checks; no statistical analysis is rerun."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def count_tsv(relative: str) -> int:
    path = ROOT / relative
    if not path.is_file():
        raise AssertionError(f"missing file: {relative}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


EXPECTED = {
    "results/figure_source_data/figure2/F2A_factor_loadings.tsv": 241,
    "results/figure_source_data/figure2/F2B_RBO_matrix.tsv": 12,
    "results/figure_source_data/figure2/F2E_factor_correspondence.tsv": 7,
    "results/summary_tables/PHASE18_HUNT_44_TEST_FINAL.tsv": 44,
    "results/summary_tables/PHASE18_UNIFIED_22_CELL_STATUS.tsv": 22,
    "results/summary_tables/PHASE18_DISEASE_FINAL_RESULTS.tsv": 40,
}

for relative, expected in EXPECTED.items():
    observed = count_tsv(relative)
    if observed != expected:
        raise AssertionError(f"{relative}: expected {expected} rows, observed {observed}")

print(f"PASS: {len(EXPECTED)} frozen-table dimensions verified")
