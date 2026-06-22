#!/usr/bin/env python3
"""
Generate all figures for the Bioinformatics submission package.

Reads bundled CSV/JSON from ../data/ (self-contained for Overleaf upload).
Outputs PNG files into this directory.

Usage:
    python figures/generate_figures.py
"""

from __future__ import annotations

import csv
import json
import os
from collections import defaultdict
from typing import Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))

PAL = {
    "error_prone": "#EF4444",
    "clean": "#0D9488",
    "muted": "#94A3B8",
    "navy": "#1E293B",
    "amber": "#F59E0B",
    "xray": "#3B82F6",
    "nmr": "#8B5CF6",
    "cryoem": "#10B981",
}

_DEPOSITION_YEARS = {
    "1ABI": 1992, "1BG0": 1998, "1D7T": 1999, "1HHZ": 2000,
    "1KO0": 2001, "1MCB": 1993, "1OF6": 2003, "1P52": 2003,
    "1UHG": 2003, "1XT7": 2004, "2AOU": 2005, "2ATS": 2005,
    "2RMI": 2007, "2W76": 2008, "2H9E": 2006, "3RIT": 2011,
}

_ERROR_TYPE_COLORS = {
    "Stereochem": "#EF4444",
    "CCD-Code": "#F59E0B",
    "Mislabel": "#3B82F6",
    "Borderline": "#94A3B8",
}


def _read_ccd_coverage() -> List[Dict]:
    rows: List[Dict] = []
    with open(os.path.join(DATA, "ccd_code_coverage_summary.csv")) as fp:
        for r in csv.DictReader(fp):
            r["error_rate_pct"] = float(r["error_rate_pct"])
            r["n_errors"] = int(r["n_errors"])
            r["residues_checked"] = int(r["residues_checked"])
            rows.append(r)
    return rows


def _read_error_table() -> List[Dict]:
    rows: List[Dict] = []
    with open(os.path.join(DATA, "error_table_verified.csv")) as fp:
        for r in csv.DictReader(fp):
            rows.append(r)
    return rows


def _read_signed_volumes() -> Dict[str, np.ndarray]:
    d_vals: List[float] = []
    l_vals: List[float] = []
    with open(os.path.join(DATA, "d_residue_verification.csv")) as fp:
        for r in csv.DictReader(fp):
            if r["has_all_atoms"] != "True":
                continue
            try:
                v = float(r["signed_volume"])
            except (TypeError, ValueError):
                continue
            if r["chirality"] == "D":
                d_vals.append(v)
            else:
                l_vals.append(v)
    return {"d": np.asarray(d_vals), "l": np.asarray(l_vals)}


def figure1_error_rate_per_ccd(rows: List[Dict]) -> str:
    rows_sorted = sorted(rows, key=lambda r: -r["error_rate_pct"])
    codes = [r["ccd_code"] for r in rows_sorted]
    rates = [r["error_rate_pct"] for r in rows_sorted]
    colors = [PAL["error_prone"] if r["n_errors"] > 0 else PAL["clean"] for r in rows_sorted]

    fig, ax = plt.subplots(figsize=(10, 4.5))
    bars = ax.bar(codes, rates, color=colors, edgecolor=PAL["navy"], linewidth=0.6)
    ax.set_xlabel("D-amino acid CCD code")
    ax.set_ylabel("Error rate (%)")
    ax.set_title(
        "Per-CCD-code D-label/L-coordinate mismatch rate\n"
        "(red = codes with errors, teal = confirmed clean at ≥91% RCSB coverage)"
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bar, r in zip(bars, rows_sorted):
        if r["n_errors"] > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.05,
                f'{r["n_errors"]}/{r["residues_checked"]}',
                ha="center", va="bottom", fontsize=8, color=PAL["navy"],
            )
    fig.tight_layout()
    out = os.path.join(HERE, "fig1_error_rate_per_ccd.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def figure2_deposition_year(rows: List[Dict]) -> str:
    counts: Dict[str, int] = defaultdict(int)
    for r in rows:
        counts[r["PDB"]] += 1

    xs, ys, labels = [], [], []
    for pdb, n in counts.items():
        year = _DEPOSITION_YEARS.get(pdb)
        if year is None:
            continue
        xs.append(year)
        ys.append(n)
        labels.append(pdb)

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.scatter(
        xs, ys, s=[80 + 40 * n for n in ys],
        color=PAL["error_prone"], alpha=0.75,
        edgecolor=PAL["navy"], linewidth=0.8,
    )
    for x, y, lbl in zip(xs, ys, labels):
        ax.annotate(lbl, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.axvspan(2006, 2008, color=PAL["muted"], alpha=0.15)
    ax.axvline(2006, linestyle="--", color=PAL["muted"], linewidth=0.8)
    ax.axvline(2008, linestyle="--", color=PAL["muted"], linewidth=0.8)
    ax.text(2007, max(ys) * 0.95, "wwPDB remediation\n2006–2008",
            ha="center", fontsize=8, color=PAL["muted"])
    ax.set_xlabel("Deposition year")
    ax.set_ylabel("D-label/L-coordinate mismatches per structure")
    ax.set_title("D-chirality errors vs deposition year")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(min(xs) - 1, max(xs) + 1)
    ax.set_ylim(0, max(ys) + 2)
    fig.tight_layout()
    out = os.path.join(HERE, "fig2_deposition_year_vs_errors.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def figure3_signed_volume(vols: Dict[str, np.ndarray]) -> str:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bins = np.linspace(-4.0, 4.0, 81)
    if vols["d"].size:
        ax.hist(vols["d"], bins=bins, color=PAL["clean"], alpha=0.75,
                label=f'D-chirality (n={vols["d"].size:,})', edgecolor="none")
    if vols["l"].size:
        ax.hist(vols["l"], bins=bins, color=PAL["error_prone"], alpha=0.85,
                label=f'L-chirality found (n={vols["l"].size})',
                edgecolor=PAL["navy"], linewidth=0.4)
    ax.axvline(0.0, color=PAL["navy"], linewidth=0.8)
    ax.set_yscale("log")
    ax.set_xlabel(r"Signed C$_\alpha$ tetrahedron volume (Å$^3$)")
    ax.set_ylabel("Residue count (log scale)")
    ax.set_title(
        "Signed-volume distribution at D-labeled C$_\\alpha$ centers"
        "\n(n = 12,573 residues; 29 with positive volume = L-coordinates)"
    )
    ax.legend(loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = os.path.join(HERE, "fig3_signed_volume_distribution.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def figure4_bland_altman() -> str:
    with open(os.path.join(DATA, "molprobity_comparison.json")) as fp:
        rows = json.load(fp)

    cf = np.array([r["cf_rama_outlier"] for r in rows], dtype=float)
    ww = np.array([r["wwpdb_rama_outlier"] for r in rows], dtype=float)
    mean = (cf + ww) / 2.0
    diff = cf - ww
    bias = float(np.mean(diff))
    sd = float(np.std(diff, ddof=1))
    loa_lo, loa_hi = bias - 1.96 * sd, bias + 1.96 * sd

    method_colors = []
    for r in rows:
        m = (r.get("method") or "").upper()
        if "NMR" in m:
            method_colors.append(PAL["nmr"])
        elif "ELECTRON" in m or "CRYO" in m:
            method_colors.append(PAL["cryoem"])
        else:
            method_colors.append(PAL["xray"])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(mean, diff, c=method_colors, s=55, alpha=0.85, edgecolor=PAL["navy"], linewidth=0.5)
    ax.axhline(bias, color=PAL["navy"], linewidth=1.2, label=f"Bias = {bias:.3f}%")
    ax.axhline(loa_lo, color=PAL["error_prone"], linestyle="--", linewidth=0.9,
               label=f"LoA = [{loa_lo:.2f}, {loa_hi:.2f}]%")
    ax.axhline(loa_hi, color=PAL["error_prone"], linestyle="--", linewidth=0.9)
    se = sd / np.sqrt(len(diff)) * 1.96
    ax.fill_between(
        [mean.min() - 0.5, mean.max() + 0.5],
        bias - se, bias + se,
        color=PAL["muted"], alpha=0.2,
    )
    ax.set_xlabel("Mean Ramachandran outlier % (ChiralFold + wwPDB) / 2")
    ax.set_ylabel("Difference (ChiralFold − wwPDB) %")
    ax.set_title(f"Bland–Altman agreement (n={len(rows)} pilot structures)")
    ax.legend(loc="upper right", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = os.path.join(HERE, "fig4_bland_altman.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def figure5_ccd_heatmap(ccd_rows: List[Dict], err_rows: List[Dict]) -> str:
    error_types = ["Stereochem", "CCD-Code", "Mislabel", "Borderline"]
    counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for r in err_rows:
        ccd = r["Residue"].split("-")[0].strip()
        if len(ccd) != 3:
            ccd = r["Residue"][:3]
        etype = r["Error_Type"]
        counts[ccd][etype] += 1

    codes_sorted = sorted(ccd_rows, key=lambda r: -r["error_rate_pct"])
    codes = [r["ccd_code"] for r in codes_sorted]
    matrix = np.array([[counts[c].get(t, 0) for t in error_types] for c in codes])
    rates = [r["error_rate_pct"] for r in codes_sorted]
    bar_colors = [PAL["error_prone"] if r["n_errors"] > 0 else PAL["clean"] for r in codes_sorted]

    fig, (ax_bar, ax_hm) = plt.subplots(
        1, 2, figsize=(12, 6.5), gridspec_kw={"width_ratios": [1.1, 1.4]}
    )
    y = np.arange(len(codes))
    ax_bar.barh(y, rates, color=bar_colors, edgecolor=PAL["navy"], linewidth=0.4)
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(codes, fontsize=9)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Error rate (%)")
    ax_bar.set_title("Per-CCD error rate")
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    im = ax_hm.imshow(matrix, aspect="auto", cmap="YlOrRd")
    ax_hm.set_xticks(range(len(error_types)))
    ax_hm.set_xticklabels(error_types, rotation=30, ha="right", fontsize=9)
    ax_hm.set_yticks(range(len(codes)))
    ax_hm.set_yticklabels(codes, fontsize=9)
    ax_hm.set_title("Error counts by type and CCD code")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = int(matrix[i, j])
            if val > 0:
                ax_hm.text(j, i, str(val), ha="center", va="center", fontsize=8, color=PAL["navy"])
    fig.colorbar(im, ax=ax_hm, fraction=0.03, pad=0.02)
    fig.tight_layout()
    out = os.path.join(HERE, "fig5_ccd_heatmap.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def figure6_ramachandran_expanded() -> str:
    rows = []
    with open(os.path.join(DATA, "ramachandran_279struct_chainfix_comparison.csv")) as fp:
        for r in csv.DictReader(fp):
            if r["pdb_id"] == "5M2K":
                continue
            rows.append(r)

    cf = np.array([float(r["chiralfold_rama_outlier_pct"]) for r in rows])
    ww = np.array([float(r["wwpdb_rama_outlier_pct"]) for r in rows])
    rho, rho_p = stats.spearmanr(cf, ww)

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.scatter(ww, cf, s=24, alpha=0.72, color=PAL["xray"], edgecolor=PAL["navy"], linewidth=0.3)
    lim = max(cf.max(), ww.max(), 1.0) * 1.08
    ax.plot([0, lim], [0, lim], "--", linewidth=0.9, color=PAL["muted"], label="y = x")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("wwPDB/MolProbity Ramachandran outlier (%)")
    ax.set_ylabel("ChiralFold Ramachandran outlier (%)")
    ax.set_title(
        f"Ramachandran agreement (n={len(rows)} structures)\n"
        f"Spearman ρ = {rho:.2f} (p = {rho_p:.2g})"
    )
    ax.legend(loc="upper left", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = os.path.join(HERE, "fig6_ramachandran_expanded.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def main() -> int:
    ccd_rows = _read_ccd_coverage()
    err_rows = _read_error_table()
    sv = _read_signed_volumes()

    outputs = [
        figure1_error_rate_per_ccd(ccd_rows),
        figure2_deposition_year(err_rows),
        figure3_signed_volume(sv),
        figure4_bland_altman(),
        figure5_ccd_heatmap(ccd_rows, err_rows),
        figure6_ramachandran_expanded(),
    ]
    for path in outputs:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
