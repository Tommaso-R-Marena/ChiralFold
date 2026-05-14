#!/usr/bin/env python3
"""
Generate the three figures used in paper/pdb_survey_note.tex.

Inputs (all already in the repo, no network required):
    results/ccd_code_coverage_summary.csv  -> per-CCD-code error rates (Fig 1)
    results/error_table_verified.csv       -> per-error-structure metadata (Fig 2)
    results/d_residue_verification.csv     -> full signed-volume distribution (Fig 3)

Outputs:
    paper/pdb_survey_figures/fig1_error_rate_per_ccd.png
    paper/pdb_survey_figures/fig2_deposition_year_vs_errors.png
    paper/pdb_survey_figures/fig3_signed_volume_distribution.png
"""

from __future__ import annotations

import csv
import os
import re
import sys
import urllib.error
import urllib.request
from collections import defaultdict
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RESULTS = os.path.join(REPO_ROOT, "results")

PAL = {
    "error_prone": "#EF4444",
    "clean": "#0D9488",
    "muted": "#94A3B8",
    "navy": "#1E293B",
    "amber": "#F59E0B",
}


def _read_ccd_coverage() -> List[Dict]:
    rows: List[Dict] = []
    with open(os.path.join(RESULTS, "ccd_code_coverage_summary.csv")) as fp:
        for r in csv.DictReader(fp):
            r["error_rate_pct"] = float(r["error_rate_pct"])
            r["n_errors"] = int(r["n_errors"])
            r["residues_checked"] = int(r["residues_checked"])
            rows.append(r)
    return rows


def _read_error_table() -> List[Dict]:
    rows: List[Dict] = []
    with open(os.path.join(RESULTS, "error_table_verified.csv")) as fp:
        for r in csv.DictReader(fp):
            rows.append(r)
    return rows


def _read_signed_volumes() -> Dict[str, np.ndarray]:
    """Return {"d": all-D-chirality vols, "l": all-L-chirality vols}."""
    d_vals: List[float] = []
    l_vals: List[float] = []
    path = os.path.join(RESULTS, "d_residue_verification.csv")
    with open(path) as fp:
        reader = csv.DictReader(fp)
        for r in reader:
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


# Deposition years are not in error_table_verified.csv; they come from
# Table 1/2 of paper/chiralfold_paper.tex. The same 16 structures.
_DEPOSITION_YEARS = {
    "1ABI": 1992, "1BG0": 1998, "1D7T": 1999, "1HHZ": 2000,
    "1KO0": 2001, "1MCB": 1993, "1OF6": 2003, "1P52": 2003,
    "1UHG": 2003, "1XT7": 2004, "2AOU": 2005, "2ATS": 2005,
    "2RMI": 2007, "2W76": 2008, "2H9E": 2006, "3RIT": 2011,
}


def figure1_error_rate_per_ccd(rows: List[Dict]) -> str:
    rows_sorted = sorted(rows, key=lambda r: -r["error_rate_pct"])
    codes = [r["ccd_code"] for r in rows_sorted]
    rates = [r["error_rate_pct"] for r in rows_sorted]
    colors = [PAL["error_prone"] if r["n_errors"] > 0 else PAL["clean"]
              for r in rows_sorted]

    fig, ax = plt.subplots(figsize=(10, 4.5))
    bars = ax.bar(codes, rates, color=colors, edgecolor=PAL["navy"], linewidth=0.6)
    ax.set_xlabel("D-amino acid CCD code")
    ax.set_ylabel("Error rate (%)")
    ax.set_title(
        "Per-CCD-code D-label/L-coordinate mismatch rate\n"
        "(red = error-prone, teal = confirmed clean at ≥91% RCSB coverage)"
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for bar, r in zip(bars, rows_sorted):
        if r["n_errors"] > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.05,
                    f'{r["n_errors"]}/{r["residues_checked"]}',
                    ha="center", va="bottom", fontsize=8, color=PAL["navy"])

    fig.tight_layout()
    out = os.path.join(HERE, "fig1_error_rate_per_ccd.png")
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def figure2_deposition_year(rows: List[Dict]) -> str:
    """Per-structure error count (size) vs deposition year (x)."""
    counts: Dict[str, int] = defaultdict(int)
    for r in rows:
        counts[r["PDB"]] += 1

    xs: List[int] = []
    ys: List[int] = []
    labels: List[str] = []
    for pdb, n in counts.items():
        year = _DEPOSITION_YEARS.get(pdb)
        if year is None:
            continue
        xs.append(year)
        ys.append(n)
        labels.append(pdb)

    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    ax.scatter(xs, ys, s=[80 + 40 * n for n in ys],
               color=PAL["error_prone"], alpha=0.75,
               edgecolor=PAL["navy"], linewidth=0.8)
    for x, y, lbl in zip(xs, ys, labels):
        ax.annotate(lbl, (x, y), xytext=(4, 4),
                    textcoords="offset points", fontsize=8)
    ax.axvline(2006, linestyle="--", color=PAL["muted"], linewidth=0.8)
    ax.axvline(2008, linestyle="--", color=PAL["muted"], linewidth=0.8)
    ax.text(2007, max(ys) * 0.95, "wwPDB remediation\n2006–2008",
            ha="center", fontsize=8, color=PAL["muted"])
    ax.set_xlabel("Deposition year")
    ax.set_ylabel("D-label/L-coordinate mismatches per structure")
    ax.set_title("D-chirality errors persist past the 2006–2008 wwPDB remediation")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(min(xs) - 1, max(xs) + 1)
    ax.set_ylim(0, max(ys) + 2)
    ax.set_yticks(range(0, max(ys) + 2))

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


def main() -> int:
    ccd_rows = _read_ccd_coverage()
    err_rows = _read_error_table()
    sv = _read_signed_volumes()

    f1 = figure1_error_rate_per_ccd(ccd_rows)
    f2 = figure2_deposition_year(err_rows)
    f3 = figure3_signed_volume(sv)
    for path in (f1, f2, f3):
        print(f"wrote {os.path.relpath(path, REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
