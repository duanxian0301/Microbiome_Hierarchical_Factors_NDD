from __future__ import annotations

import bisect
import csv
import gzip
import json
import math
from collections import defaultdict
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ["MHFNDD_PROJECT_ROOT"]).resolve()
BASE = ROOT / "FULL_SCALE_GFA_FINAL_EVALUATION"
P11 = ROOT / "PHASE11_DISEASE_BIOLOGICAL_INTEGRATION"
P11X = P11 / "EXPLORATORY_STANDARD_GWAS_DOWNSTREAM"
OUT = ROOT / "PHASE13_SHARED_STRUCTURE_AND_DISEASE_LOCI" / "B_FOUR_DISEASE_LOCI"
FIG = OUT / "PHASE13B_FIGURE_CANDIDATE_LOCI"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
REF = ROOT / "GLOBAL_LOCAL_HETEROGENEITY_CELL_CONTEXT" / "CELL_ANNOTATION_C1_EXECUTION" / "reference_download" / "GRCh38" / "plink_files"
DISEASES = {"AD", "SCZ", "CD", "UC"}


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def p_from_z(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0)) if math.isfinite(float(z)) else float("nan")


def bh(pvals: list[float]) -> list[float]:
    n = len(pvals)
    order = np.argsort(np.asarray(pvals, float))
    out = np.ones(n, float)
    running = 1.0
    for rank0 in range(n - 1, -1, -1):
        idx = int(order[rank0])
        running = min(running, pvals[idx] * n / (rank0 + 1))
        out[idx] = min(1.0, running)
    return out.tolist()


def norm_key(key: str) -> tuple[str, int, str, str]:
    c, p, a1, a2 = key.split(":")
    a, b = sorted((a1, a2))
    return c, int(p), a, b


def load_regions() -> list[dict]:
    src = P11X / "04_DISEASE_REGION_LOOKUP" / "DISEASE_REGION_M1_LOOKUP.tsv"
    rows = []
    with src.open(encoding="utf-8") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r["DISEASE"] not in DISEASES:
                continue
            x = dict(r)
            for k in ("CHR", "START", "END", "DISEASE_LEAD_BP", "N_CLUMP_LEADS_MERGED"):
                x[k] = int(float(x[k]))
            x["DISEASE_LEAD_P"] = float(x["DISEASE_LEAD_P"])
            rows.append(x)
    counts = defaultdict(int)
    for r in rows:
        counts[r["DISEASE"]] += 1
    if dict(counts) != {"SCZ": 130, "AD": 41, "CD": 77, "UC": 48}:
        raise ValueError(counts)
    return rows


def load_bim() -> tuple[dict[str, tuple[str, int, str, str, tuple[str, int]]], dict[tuple[str, int, str, str], tuple[str, int]]]:
    rs = {}
    key = {}
    for chrom in range(1, 23):
        with (REF / f"1000G.EUR.hg38.{chrom}.bim").open() as handle:
            for i, line in enumerate(handle):
                c, rid, _, bp, a1, a2 = line.split()
                a, b = sorted((a1, a2))
                bed_index = (c, i)
                rs[rid] = (c, int(bp), a, b, bed_index)
                key[(c, int(bp), a, b)] = bed_index
    return rs, key


def build_region_indices(keys: np.ndarray, regions: list[dict], lead_tuples: dict[tuple[str, int, str, str], list[int]]) -> tuple[list[np.ndarray], dict[int, int]]:
    by_chr = defaultdict(list)
    for i, r in enumerate(regions):
        by_chr[int(r["CHR"])].append((int(r["START"]), int(r["END"]), i))
    for c in by_chr:
        by_chr[c].sort()
    out = [[] for _ in regions]
    lead_index: dict[int, int] = {}
    current_chr = None
    ptr = 0
    active: list[tuple[int, int]] = []
    intervals = []
    for idx, raw in enumerate(keys):
        s = str(raw)
        c_s, p_s, a1, a2 = s.split(":")
        c = int(c_s)
        pos = int(p_s)
        if c != current_chr:
            current_chr = c
            ptr = 0
            active = []
            intervals = by_chr.get(c, [])
        while ptr < len(intervals) and intervals[ptr][0] <= pos:
            active.append((intervals[ptr][1], intervals[ptr][2]))
            ptr += 1
        active = [(end, ri) for end, ri in active if end >= pos]
        for _, ri in active:
            out[ri].append(idx)
        a, b = sorted((a1, a2))
        nk = (c_s, pos, a, b)
        for ri in lead_tuples.get(nk, []):
            lead_index[ri] = idx
    return [np.asarray(x, dtype=np.int64) for x in out], lead_index


def genotype_reader():
    nind = 489
    nbyte = (nind + 3) // 4
    lut = np.empty((256, 4), dtype=np.float32)
    code = [0.0, np.nan, 1.0, 2.0]
    for byte in range(256):
        for j in range(4):
            lut[byte, j] = code[(byte >> (2 * j)) & 3]
    handles = {}
    for chrom in range(1, 23):
        handle = (REF / f"1000G.EUR.hg38.{chrom}.bed").open("rb")
        if handle.read(3) != b"\x6c\x1b\x01":
            raise ValueError(f"Unexpected BED header for chromosome {chrom}")
        handles[str(chrom)] = handle
    cache: dict[tuple[str, int], np.ndarray] = {}

    def geno(i: tuple[str, int] | None) -> np.ndarray | None:
        if i is None:
            return None
        if i not in cache:
            chrom, row = i
            handle = handles[chrom]
            handle.seek(3 + row * nbyte)
            cache[i] = lut[np.frombuffer(handle.read(nbyte), dtype=np.uint8)].reshape(-1)[:nind].copy()
        return cache[i]

    def r2(i: tuple[str, int] | None, j: tuple[str, int] | None) -> float | None:
        if i is not None and j is not None and i[0] != j[0]:
            return 0.0
        a = geno(i)
        b = geno(j)
        if a is None or b is None:
            return None
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 30 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
            return None
        return float(np.corrcoef(a[ok], b[ok])[0, 1] ** 2)

    return r2, handles


def peak_context(candidate_key: str, region: dict, beds: dict[str, Path]) -> list[dict]:
    c, pos, _, _ = norm_key(candidate_key)
    rows = []
    for cell, path in beds.items():
        overlap = False
        region_peaks = 0
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                v = line.rstrip("\n").split("\t")
                cc = v[0].removeprefix("chr")
                if cc != c:
                    continue
                start0, end = int(v[1]), int(v[2])
                if start0 < pos <= end:
                    overlap = True
                if end >= int(region["START"]) and start0 + 1 <= int(region["END"]):
                    region_peaks += 1
        rows.append({"CELL": cell, "CANDIDATE_SNP_IN_REGULATORY_PEAK": overlap, "REGION_OVERLAPPING_PEAK_N": region_peaks})
    return rows


def make_svg(path: Path, candidates: list[dict]) -> None:
    top = candidates[:12]
    w, h = 1100, 520
    left, top_y, pw, ph = 125, 55, 920, 350
    vmax = max(float(x["PRIORITY_SCORE"]) for x in top) if top else 1
    bw = pw / max(1, len(top)) * 0.65
    colors = {"AD": "#4472c4", "SCZ": "#7030a0", "CD": "#70ad47", "UC": "#ed7d31"}
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">', '<rect width="100%" height="100%" fill="white"/>', '<text x="24" y="28" font-family="Arial" font-size="18" font-weight="bold">Phase13B prioritized disease-region and core-species intersections</text>']
    for i, x in enumerate(top):
        score = float(x["PRIORITY_SCORE"])
        xx = left + (i + 0.5) * pw / len(top) - bw / 2
        hh = ph * score / vmax
        yy = top_y + ph - hh
        parts.append(f'<rect x="{xx:.1f}" y="{yy:.1f}" width="{bw:.1f}" height="{hh:.1f}" fill="{colors.get(x["DISEASE"],"#888")}"/>')
        label = f'{x["DISEASE"]} chr{x["CHR"]}:{int(x["START"])//1000000}Mb'
        parts.append(f'<text x="{xx+bw/2:.1f}" y="{top_y+ph+18}" text-anchor="end" transform="rotate(-55 {xx+bw/2:.1f} {top_y+ph+18})" font-family="Arial" font-size="10">{label}</text>')
        parts.append(f'<text x="{xx+bw/2:.1f}" y="{yy-5:.1f}" text-anchor="middle" font-family="Arial" font-size="10">{score:.2f}</text>')
    parts.append('<text x="28" y="250" transform="rotate(-90 28 250)" font-family="Arial" font-size="12">Priority score</text>')
    parts.append('<text x="125" y="505" font-family="Arial" font-size="11" fill="#555">Scores rank follow-up value; they are not calibrated probabilities or evidence of colocalization.</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    regions = load_regions()
    rs_map, bim_key = load_bim()
    source_manifest = pd.read_csv(BASE / "inputs_241_manifest.tsv", sep="\t")
    h2 = pd.read_csv(BASE / "BIDIRECTIONAL_SPECIES_LOADING_H2.tsv", sep="\t")
    load = pd.read_csv(BASE / "cross_cohort_validation" / "aligned_species_loadings.tsv", sep="\t")
    species = h2[(h2["Swedish_F1_top40"] == True) | (h2["Swedish_F2_top40"] == True)].copy()
    species["PRIMARY_TOP20_UNION"] = (species["Swedish_F1_top20"] == True) | (species["Swedish_F2_top20"] == True)
    species = species.sort_values(["PRIMARY_TOP20_UNION", "trait_order"], ascending=[False, True])
    frozen_rows = []
    for _, r in species.iterrows():
        frozen_rows.append({
            "TRAIT_ORDER": int(r["trait_order"]), "SPECIES": r["species"], "SWEDISH_ACCESSION": r["swedish_accession"], "HUNT_ACCESSION": r["hunt_accession"],
            "PRIMARY_TOP20_UNION": bool(r["PRIMARY_TOP20_UNION"]), "EXPANDED_TOP40_UNION": True,
            "F1_LOADING": r["Swedish_F1_loading_raw"], "F1_ABS_RANK": int(r["Swedish_F1_abs_rank"]), "F1_SIGN": "positive" if r["Swedish_F1_loading_raw"] > 0 else "negative",
            "F2_LOADING": r["Swedish_F2_loading_raw"], "F2_ABS_RANK": int(r["Swedish_F2_abs_rank"]), "F2_SIGN": "positive" if r["Swedish_F2_loading_raw"] > 0 else "negative",
            "SWEDISH_H2": r["Swedish_h2"], "SWEDISH_H2_SE": r["Swedish_h2_SE"], "SWEDISH_H2_Z": r["Swedish_h2_Z"],
            "ORIGINAL_GWAS_N": 16017, "QC_AVAILABLE": "beta,se,z; no complete per-variant INFO/MAF in dense NPZ",
        })
    write_tsv(OUT / "PHASE13B_FROZEN_SPECIES_LIST.tsv", frozen_rows)

    disease_sources = {
        "AD": ("GCST90704648 clinical no-biobank European", "GRCh38", "41"),
        "SCZ": ("PGC3 EUR wave 3", "GRCh38 harmonized", "130"),
        "CD": ("de Lange 2017 GCST004132", "GRCh37 lifted to GRCh38", "77"),
        "UC": ("de Lange 2017 GCST004133", "GRCh37 lifted to GRCh38", "48"),
    }
    region_manifest = []
    for r in regions:
        src, build, n = disease_sources[r["DISEASE"]]
        region_manifest.append({
            "REGION_ID": r["REGION_ID"], "DISEASE": r["DISEASE"], "CHR": r["CHR"], "START": r["START"], "END": r["END"],
            "DISEASE_LEAD_SNP": r["DISEASE_LEAD_SNP"], "DISEASE_LEAD_BP": r["DISEASE_LEAD_BP"], "DISEASE_LEAD_P": r["DISEASE_LEAD_P"],
            "N_CLUMP_LEADS_MERGED": r["N_CLUMP_LEADS_MERGED"], "ALL_CLUMP_LEADS": r["ALL_CLUMP_LEADS"], "NEAREST_GENE": r["NEAREST_GENE"],
            "MHC_REGION": r["MHC_REGION"], "AD_APOE_REGION": r["AD_APOE_REGION"], "SOURCE": src, "SOURCE_BUILD": build,
            "REGION_DEFINITION": "disease-only HM3 clump P1=5e-8 r2=0.1 1Mb; lead +/-250kb; merge overlap/gap<=100kb", "FROZEN_DISEASE_REGION_COUNT": n,
        })
    write_tsv(OUT / "PHASE13B_DISEASE_REGION_MANIFEST.tsv", region_manifest)

    first_npz = Path(source_manifest.iloc[0]["swedish_matrix"])
    keys = np.load(first_npz, allow_pickle=False)["snp"]
    lead_tuples: dict[tuple[str, int, str, str], list[int]] = defaultdict(list)
    lead_bim_row = {}
    for ri, r in enumerate(regions):
        info = rs_map.get(r["DISEASE_LEAD_SNP"])
        if info:
            c, bp, a, b, bi = info
            lead_tuples[(str(r["CHR"]), bp, a, b)].append(ri)
            lead_bim_row[ri] = bi
    region_indices, lead_indices = build_region_indices(keys, regions, lead_tuples)
    if min(len(x) for x in region_indices) == 0:
        raise ValueError("At least one frozen region has no dense species SNP")
    checkpoint = {"status": "region_index_ready", "regions": len(regions), "species_total": len(species), "lead_indices": len(lead_indices), "completed_species": []}
    (OUT / "checkpoint.json").write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")

    accession_to_manifest = source_manifest.set_index("swedish_accession")
    all_rows: list[dict] = []
    for si, (_, sp) in enumerate(species.iterrows(), 1):
        accession = sp["swedish_accession"]
        npz_path = Path(accession_to_manifest.loc[accession, "swedish_matrix"])
        data = np.load(npz_path, allow_pickle=False)
        z = data["z"]
        beta = data["beta"]
        se = data["se"]
        for ri, r in enumerate(regions):
            idxs = region_indices[ri]
            zz = z[idxs]
            finite = np.isfinite(zz)
            if not finite.any():
                continue
            local = int(np.nanargmax(np.abs(zz)))
            min_idx = int(idxs[local])
            min_z = float(z[min_idx])
            lead_idx = lead_indices.get(ri)
            lead_z = float(z[lead_idx]) if lead_idx is not None else float("nan")
            all_rows.append({
                "REGION_ID": r["REGION_ID"], "DISEASE": r["DISEASE"], "CHR": r["CHR"], "START": r["START"], "END": r["END"],
                "DISEASE_LEAD_SNP": r["DISEASE_LEAD_SNP"], "DISEASE_LEAD_BP": r["DISEASE_LEAD_BP"], "DISEASE_LEAD_P": r["DISEASE_LEAD_P"], "AD_APOE_REGION": r["AD_APOE_REGION"],
                "SPECIES": sp["species"], "SWEDISH_ACCESSION": accession, "PRIMARY_TOP20_UNION": bool(sp["PRIMARY_TOP20_UNION"]),
                "F1_LOADING": sp["Swedish_F1_loading_raw"], "F1_ABS_RANK": int(sp["Swedish_F1_abs_rank"]), "F2_LOADING": sp["Swedish_F2_loading_raw"], "F2_ABS_RANK": int(sp["Swedish_F2_abs_rank"]),
                "LEAD_SNP_DENSE_KEY": str(keys[lead_idx]) if lead_idx is not None else "NA", "LEAD_SPECIES_BETA": float(beta[lead_idx]) if lead_idx is not None else "NA", "LEAD_SPECIES_SE": float(se[lead_idx]) if lead_idx is not None else "NA",
                "LEAD_SPECIES_Z": lead_z if lead_idx is not None else "NA", "LEAD_SPECIES_P": p_from_z(lead_z) if lead_idx is not None else "NA", "LEAD_SPECIES_EFFECT_DIRECTION": "positive" if lead_z > 0 else ("negative" if lead_z < 0 else "NA"),
                "REGIONAL_MIN_P": p_from_z(min_z), "REGIONAL_MIN_Z": min_z, "REGIONAL_MIN_BETA": float(beta[min_idx]), "REGIONAL_MIN_SE": float(se[min_idx]), "REGIONAL_MIN_KEY": str(keys[min_idx]),
                "REGION_DENSE_SNP_N": int(finite.sum()), "LEAD_TO_REGIONAL_MIN_R2": "PENDING", "REGION_REPORTED_IN_SWEDISH_OR_HUNT": "PENDING_NOVELTY_CROSSWALK",
                "INTERPRETATION": "regional minimum is a screening statistic; LD and multiplicity must be considered",
            })
        checkpoint["completed_species"].append(str(sp["species"]))
        checkpoint["status"] = f"species_{si}_of_{len(species)}"
        if si % 5 == 0 or si == len(species):
            (OUT / "checkpoint.json").write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")
    del z, beta, se

    # Empirical EUR LD between each disease lead and each species regional minimum.
    candidate_tuples = {norm_key(r["REGIONAL_MIN_KEY"]) for r in all_rows}
    candidate_bim = {k: bim_key.get(k) for k in candidate_tuples}
    r2fun, bed_handles = genotype_reader()
    region_id_to_i = {x["REGION_ID"]: i for i, x in enumerate(regions)}
    for r in all_rows:
        ri = region_id_to_i[r["REGION_ID"]]
        val = r2fun(lead_bim_row.get(ri), candidate_bim.get(norm_key(r["REGIONAL_MIN_KEY"])))
        r["LEAD_TO_REGIONAL_MIN_R2"] = round(val, 6) if val is not None else "NA"
    for bed_handle in bed_handles.values():
        bed_handle.close()

    df = pd.DataFrame(all_rows)
    primary = df[df["PRIMARY_TOP20_UNION"] == True].copy()
    agg_rows = []
    for ri, region in enumerate(regions):
        q = primary[primary["REGION_ID"] == region["REGION_ID"]].copy()
        q["R2NUM"] = pd.to_numeric(q["LEAD_TO_REGIONAL_MIN_R2"], errors="coerce")
        best = q.sort_values(["REGIONAL_MIN_P", "R2NUM"], ascending=[True, False]).iloc[0]
        leadp = pd.to_numeric(q["LEAD_SPECIES_P"], errors="coerce")
        n_nom = int((q["REGIONAL_MIN_P"] < 0.05).sum())
        n_1e4 = int((q["REGIONAL_MIN_P"] < 1e-4).sum())
        n_ld = int(((q["REGIONAL_MIN_P"] < 1e-4) & (q["R2NUM"] >= 0.1)).sum())
        n_lead_nom = int((leadp < 0.05).sum())
        best_lead_p = float(leadp.min()) if leadp.notna().any() else 1.0
        max_r2 = float(q["R2NUM"].max()) if q["R2NUM"].notna().any() else 0.0
        score = min(12, -math.log10(float(region["DISEASE_LEAD_P"]))) / 4 + min(8, -math.log10(float(best["REGIONAL_MIN_P"]))) + 2 * max_r2 + min(4, n_ld) + min(3, n_lead_nom / 3)
        if n_ld >= 2 and float(best["REGIONAL_MIN_P"]) < 1e-4:
            tier = "A"
        elif (float(best["REGIONAL_MIN_P"]) < 1e-5 and max_r2 >= 0.1) or best_lead_p < 1e-4:
            tier = "B"
        elif n_nom >= 3 or float(best["REGIONAL_MIN_P"]) < 1e-4:
            tier = "C"
        else:
            tier = "NOT_PRIORITIZED"
        agg_rows.append({
            "REGION_ID": region["REGION_ID"], "DISEASE": region["DISEASE"], "CHR": region["CHR"], "START": region["START"], "END": region["END"],
            "DISEASE_LEAD_SNP": region["DISEASE_LEAD_SNP"], "DISEASE_LEAD_P": region["DISEASE_LEAD_P"], "NEAREST_GENE": region["NEAREST_GENE"], "AD_APOE_REGION": region["AD_APOE_REGION"],
            "BEST_SPECIES": best["SPECIES"], "BEST_SPECIES_REGIONAL_MIN_P": best["REGIONAL_MIN_P"], "BEST_SPECIES_REGIONAL_MIN_KEY": best["REGIONAL_MIN_KEY"], "BEST_SPECIES_REGIONAL_MIN_Z": best["REGIONAL_MIN_Z"],
            "BEST_LEAD_SPECIES_P": best_lead_p, "MAX_LEAD_TO_SPECIES_CANDIDATE_R2": max_r2, "CORE_SPECIES_NOMINAL_REGION_N": n_nom, "CORE_SPECIES_P_LT_1E4_N": n_1e4,
            "LD_COMPATIBLE_P_LT_1E4_SPECIES_N": n_ld, "CORE_SPECIES_LEAD_P_LT_0P05_N": n_lead_nom, "TIER": tier, "PRIORITY_SCORE": score,
            "STATISTICAL_SCOPE": "screening tier; correlated species and searched SNPs; not colocalization",
        })
    # Matched search-volume background: empirical tail within region-SNP-count quartile.
    adf = pd.DataFrame(agg_rows)
    region_n = primary.groupby("REGION_ID")["REGION_DENSE_SNP_N"].first()
    adf["REGION_DENSE_SNP_N"] = adf["REGION_ID"].map(region_n)
    adf["SNP_COUNT_QUARTILE"] = pd.qcut(adf["REGION_DENSE_SNP_N"], 4, labels=False, duplicates="drop")
    tailp = []
    for _, x in adf.iterrows():
        bg = adf[adf["SNP_COUNT_QUARTILE"] == x["SNP_COUNT_QUARTILE"]]
        tailp.append((1 + int((bg["LD_COMPATIBLE_P_LT_1E4_SPECIES_N"] >= x["LD_COMPATIBLE_P_LT_1E4_SPECIES_N"]).sum())) / (1 + len(bg)))
    adf["MATCHED_SEARCH_VOLUME_EMPIRICAL_P"] = tailp
    adf["MATCHED_SEARCH_VOLUME_BH_Q"] = bh(tailp)
    agg_rows = adf.to_dict("records")
    # Append region summary fields to every species row.
    amap = {r["REGION_ID"]: r for r in agg_rows}
    for r in all_rows:
        a = amap[r["REGION_ID"]]
        r["CORE_SPECIES_NOMINAL_REGION_N"] = a["CORE_SPECIES_NOMINAL_REGION_N"]
        r["CORE_SPECIES_P_LT_1E4_N"] = a["CORE_SPECIES_P_LT_1E4_N"]
        r["LD_COMPATIBLE_P_LT_1E4_SPECIES_N"] = a["LD_COMPATIBLE_P_LT_1E4_SPECIES_N"]
        r["MATCHED_SEARCH_VOLUME_EMPIRICAL_P"] = a["MATCHED_SEARCH_VOLUME_EMPIRICAL_P"]
        r["MATCHED_SEARCH_VOLUME_BH_Q"] = a["MATCHED_SEARCH_VOLUME_BH_Q"]
    write_tsv(OUT / "PHASE13B_ALL_LOCUS_SPECIES_RESULTS.tsv", all_rows)

    # Keep the candidate table focused on stronger screening tiers. Tier C and
    # non-prioritized regions remain fully represented in the all-results file.
    candidates = [r for r in agg_rows if r["TIER"] in {"A", "B"}]
    candidates.sort(key=lambda x: (-float(x["PRIORITY_SCORE"]), x["DISEASE"], int(x["CHR"]), int(x["START"])))
    if not candidates:
        candidates = sorted(agg_rows, key=lambda x: -float(x["PRIORITY_SCORE"]))[:10]
        for r in candidates:
            r["TIER"] = "C_DESCRIPTIVE_FALLBACK"

    beds = {
        "BNGA": BASE / "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA" / "annotations" / "brain" / "CATlas2023_BNGA_GRCh38_autosomes.bed.gz",
        "SNCG_1": BASE / "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA" / "annotations" / "brain" / "CATlas2023_SNCG_1_GRCh38_autosomes.bed.gz",
        "SMC": BASE / "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA" / "annotations" / "brain" / "CATlas2023_SMC_GRCh38_autosomes.bed.gz",
        "Macrophage": BASE / "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA" / "annotations" / "gut" / "Zhao2026_Macrophage_GRCh38_autosomes.bed.gz",
        "NaiveB": BASE / "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA" / "annotations" / "gut" / "Zhao2026_NaiveB_GRCh38_autosomes.bed.gz",
        "Early_Colonocyte": BASE / "PHASE9_SYSTEMATIC_CELL_REGULATORY_GFA" / "annotations" / "gut" / "Zhao2026_Early_Colonocyte_GRCh38_autosomes.bed.gz",
    }
    phase12 = {"SNCG_1": "SCZ positive Phase12A coefficient; SCZ-SNCG_1 Swedish GLOBAL F2 Top5 negative Phase12B coefficient", "SMC": "UC positive Phase12A coefficient; CD/UC joint-model sensitivity", "Macrophage": "CD nominal Phase12A coefficient", "Early_Colonocyte": "UC nominal Phase12A coefficient", "BNGA": "no frozen disease coefficient highlight", "NaiveB": "CD nominal Phase12A coefficient"}
    cell_rows = []
    for c in candidates[:30]:
        for x in peak_context(c["BEST_SPECIES_REGIONAL_MIN_KEY"], c, beds):
            cell_rows.append({
                "REGION_ID": c["REGION_ID"], "DISEASE": c["DISEASE"], "BEST_SPECIES": c["BEST_SPECIES"], "CANDIDATE_KEY": c["BEST_SPECIES_REGIONAL_MIN_KEY"],
                "CELL": x["CELL"], "CANDIDATE_SNP_IN_REGULATORY_PEAK": x["CANDIDATE_SNP_IN_REGULATORY_PEAK"], "REGION_OVERLAPPING_PEAK_N": x["REGION_OVERLAPPING_PEAK_N"],
                "PHASE12_CONNECTION": phase12[x["CELL"]], "EQTL_CAQTl_STATUS": "NOT_TESTED_NO_FROZEN_LOCAL_QTL_RESOURCE", "MECHANISTIC_INTERPRETATION": "peak overlap is positional context only",
            })
    write_tsv(OUT / "PHASE13B_CELL_REGULATORY_CONTEXT.tsv", cell_rows)
    cell_count = defaultdict(int)
    for r in cell_rows:
        if str(r["CANDIDATE_SNP_IN_REGULATORY_PEAK"]).lower() == "true":
            cell_count[r["REGION_ID"]] += 1
    for c in candidates:
        c["CANDIDATE_SNP_REGULATORY_CELL_N"] = cell_count.get(c["REGION_ID"], 0)
        if c["CANDIDATE_SNP_REGULATORY_CELL_N"]:
            c["PRIORITY_SCORE"] = float(c["PRIORITY_SCORE"]) + 0.25 * c["CANDIDATE_SNP_REGULATORY_CELL_N"]
    candidates.sort(key=lambda x: -float(x["PRIORITY_SCORE"]))
    write_tsv(OUT / "PHASE13B_PRIORITIZED_CANDIDATES.tsv", candidates)

    feas_rows = []
    for c in candidates[:10]:
        fut2 = c["DISEASE"] == "CD" and int(c["CHR"]) == 19 and int(c["START"]) <= 48703417 <= int(c["END"])
        feas_rows.append({
            "REGION_ID": c["REGION_ID"], "DISEASE": c["DISEASE"], "BEST_SPECIES": c["BEST_SPECIES"], "CANDIDATE_KEY": c["BEST_SPECIES_REGIONAL_MIN_KEY"],
            "SPECIES_BETA_SE": "AVAILABLE", "SPECIES_MAF": "REFERENCE_ONLY_NOT_ORIGINAL", "SPECIES_SNP_SPECIFIC_N": "NOT_AVAILABLE", "SPECIES_INFO": "NOT_AVAILABLE",
            "DISEASE_BETA_SE": "AVAILABLE_IN_FROZEN_HARMONIZED_RESOURCE_FOR_HM3; FULL_REGION_REQUIRES_PREPARATION", "DISEASE_MAF_N": "DISEASE_DEPENDENT_PARTIAL", "LD_REFERENCE": "1000G_EUR_489_AVAILABLE",
            "FORMAL_COLOCALIZATION_STATUS": "NOT_READY", "REASON": "original species MAF/INFO and SNP-specific N absent; full disease-region input not uniformly prepared",
            "FUT2_PHASE11E_BOUNDARY": "NO_ROBUST_LOCAL_SIGNAL; do not rerun or relabel" if fut2 else "NOT_FUT2_PHASE11E_REGION",
        })
    write_tsv(OUT / "PHASE13B_COLOCALIZATION_FEASIBILITY.tsv", feas_rows)

    published = pd.read_csv(BASE / "PHASE2_FEASIBILITY" / "PHASE2_PUBLISHED_LOCUS_SPECIES_MASTER.tsv", sep="\t")
    novelty_rows = []
    for c in candidates:
        local = []
        for _, p in published.iterrows():
            rid = str(p["lead_SNP"])
            if rid in rs_map:
                published_chr, bp, _, _, _ = rs_map[rid]
                if str(c["CHR"]) == published_chr and int(c["START"]) <= bp <= int(c["END"]):
                    local.append(f'{p["locus"]}:{rid}:{p["reported_species_or_phenotype"]}')
        fut2 = c["DISEASE"] == "CD" and int(c["CHR"]) == 19 and int(c["START"]) <= 48703417 <= int(c["END"])
        if fut2:
            cls = "1_KNOWN_HOST_SPECIES_ASSOCIATION_AND_KNOWN_DISEASE_REGION"
            note = "FUT2 host-microbiome biology is published; Phase11E remains NO_ROBUST_LOCAL_SIGNAL"
        elif local:
            cls = "4_POTENTIAL_ADDED_MULTISPECIES_ORGANIZATION_AT_KNOWN_HOST_MICROBIOME_LOCUS"
            note = "published locus overlaps region; exact candidate-species novelty requires article-level verification"
        else:
            cls = "2_KNOWN_DISEASE_ASSOCIATION_POSSIBLE_NEW_SPECIES_LINK_NOT_YET_VERIFIED"
            note = "disease region is established by construction; microbiome novelty remains unverified"
        novelty_rows.append({"REGION_ID": c["REGION_ID"], "DISEASE": c["DISEASE"], "BEST_SPECIES": c["BEST_SPECIES"], "NOVELTY_CLASS": cls, "LOCAL_PUBLISHED_HOST_MICROBIOME_MATCHES": ";".join(local) if local else "NONE_IN_FROZEN_MASTER", "NOVELTY_NOTE": note, "CLAIM_BOUNDARY": "do not call novel without external literature verification"})
    write_tsv(OUT / "PHASE13B_NOVELTY_CROSSWALK.tsv", novelty_rows)
    make_svg(FIG / "PHASE13B_FIGURE_CANDIDATE_LOCI.svg", candidates)
    checkpoint.update({"status": "complete", "all_rows": len(all_rows), "candidate_rows": len(candidates), "cell_rows": len(cell_rows), "coloc_rows": len(feas_rows)})
    (OUT / "checkpoint.json").write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
