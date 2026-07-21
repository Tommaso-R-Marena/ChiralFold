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
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

try:
    from adjustText import adjust_text

    _HAS_ADJUST_TEXT = True
except ImportError:
    _HAS_ADJUST_TEXT = False

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))

PAL = {
    "error_prone": "#EF4444",
    "clean": "#0D9488",
    "muted": "#64748B",
    "navy": "#1E293B",
    "amber": "#F59E0B",
    "xray": "#3B82F6",
    "nmr": "#0E7490",
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
    "Borderline": "#64748B",
}

_DPI = 300


def _savefig(fig: plt.Figure, out: str) -> None:
    fig.savefig(out, dpi=_DPI, bbox_inches="tight")


def _style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _wilson_ci(k: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """Wilson score 95% CI for a binomial proportion (returns fraction, not %)."""
    if n <= 0:
        return 0.0, 0.0
    p = k / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p + z2 / (2.0 * n)) / denom
    margin = z * np.sqrt((p * (1.0 - p) + z2 / (4.0 * n)) / n) / denom
    return max(0.0, center - margin), min(1.0, center + margin)


def _aggregate_error_rate_pct(rows: List[Dict]) -> float:
    total_errors = sum(r["n_errors"] for r in rows)
    total_residues = sum(r["residues_checked"] for r in rows)
    if total_residues == 0:
        return 0.0
    return 100.0 * total_errors / total_residues


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

    yerr_lo: List[float] = []
    yerr_hi: List[float] = []
    for r, rate in zip(rows_sorted, rates):
        ci_lo, ci_hi = _wilson_ci(r["n_errors"], r["residues_checked"])
        yerr_lo.append(max(0.0, rate - 100.0 * ci_lo))
        yerr_hi.append(max(0.0, 100.0 * ci_hi - rate))

    mean_rate = _aggregate_error_rate_pct(rows)

    fig, ax = plt.subplots(figsize=(10, 4.5))
    x = np.arange(len(codes))
    bars = ax.bar(
        x,
        rates,
        color=colors,
        edgecolor=PAL["navy"],
        linewidth=0.6,
        yerr=[yerr_lo, yerr_hi],
        capsize=3,
        error_kw={"elinewidth": 0.9, "ecolor": PAL["navy"], "capthick": 0.9},
    )
    ax.axhline(mean_rate, linestyle="--", color=PAL["navy"], linewidth=1.0, alpha=0.75)
    ax.set_xticks(x)
    ax.set_xticklabels(codes, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("D-amino acid CCD code")
    ax.set_ylabel("Error rate (%)")
    ymax = max(r + hi for r, hi in zip(rates, yerr_hi))
    # Headroom so labels sit above CIs, not on bars; clean (0-error) codes unlabeled.
    ax.set_ylim(0, max(ymax * 1.32, mean_rate * 2.0, 0.9))
    for idx, (bar, r) in enumerate(zip(bars, rows_sorted)):
        if r["n_errors"] <= 0:
            continue
        label = f'{r["n_errors"]}/{r["residues_checked"]}'
        y_text = bar.get_height() + yerr_hi[idx] + 0.12
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            y_text,
            label,
            ha="center",
            va="bottom",
            fontsize=7,
            color=PAL["navy"],
            clip_on=False,
            zorder=5,
            bbox={
                "boxstyle": "round,pad=0.18",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.95,
            },
        )

    _style_axes(ax)
    fig.tight_layout()
    out = os.path.join(HERE, "fig1_error_rate_per_ccd.png")
    _savefig(fig, out)
    plt.close(fig)
    return out


def _dominant_error_type(err_rows: List[Dict], pdb: str) -> str:
    types = [r.get("Error_Type", "") for r in err_rows if r["PDB"] == pdb]
    if not types:
        return "Stereochem"
    return Counter(types).most_common(1)[0][0]


def _manual_label_offsets(n: int) -> List[Tuple[float, float]]:
    """Staggered offsets (points) to reduce label overlap when adjustText is absent."""
    base = [
        (5, 5), (8, -6), (-10, 7), (12, 2), (-6, -8),
        (3, 10), (-12, -3), (10, -9), (-4, 11), (6, -11),
        (-9, 4), (11, 8), (-7, -10), (4, -7), (9, 12),
        (-11, 6),
    ]
    return [base[i % len(base)] for i in range(n)]


def _mismatch_counts_from_verification() -> Dict[str, int]:
    """Count positive signed volumes per PDB (true mismatch multiplicity)."""
    path = os.path.join(DATA, "d_residue_verification.csv")
    counts: Dict[str, int] = defaultdict(int)
    if not os.path.isfile(path):
        return counts
    with open(path) as fp:
        for r in csv.DictReader(fp):
            try:
                if float(r.get("signed_volume") or 0) > 0:
                    counts[r["pdb_id"].upper()] += 1
            except (TypeError, ValueError):
                continue
    return counts


def figure2_deposition_year(rows: List[Dict]) -> str:
    counts = _mismatch_counts_from_verification()
    if not counts:
        for r in rows:
            counts[r["PDB"]] += 1

    xs, ys, labels, point_colors = [], [], [], []
    for pdb, n in sorted(counts.items()):
        year = _DEPOSITION_YEARS.get(pdb)
        if year is None:
            continue
        xs.append(year)
        ys.append(n)
        labels.append(pdb)
        etype = _dominant_error_type(rows, pdb)
        point_colors.append(_ERROR_TYPE_COLORS.get(etype, PAL["muted"]))

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    sizes = [55 + 28 * n for n in ys]
    ax.scatter(
        xs, ys, s=sizes, c=point_colors, alpha=0.88,
        edgecolor=PAL["navy"], linewidth=0.8, zorder=3,
    )

    # Label ONLY multiplicity outliers (n>1), with long leaders clear of markers.
    highlight = [(x, y, lbl) for x, y, lbl in zip(xs, ys, labels) if y > 1]
    # Offset in data coordinates so labels stay far from large bubbles.
    manual_xy = {
        "1OF6": (1999.0, 9.35),   # far left-up, clear of the n=8 bubble
        "2ATS": (2001.2, 5.55),   # left of n=3
        "3RIT": (2013.2, 6.85),   # right-up of n=5
    }
    for x, y, lbl in highlight:
        tx, ty = manual_xy.get(lbl, (x + 1.2, y + 1.5))
        ax.annotate(
            lbl,
            (x, y),
            xytext=(tx, ty),
            textcoords="data",
            fontsize=8,
            fontweight="bold",
            color=PAL["navy"],
            ha="center",
            va="center",
            arrowprops={
                "arrowstyle": "-|>",
                "color": PAL["navy"],
                "lw": 0.75,
                "shrinkA": 3,
                "shrinkB": 14,
                "connectionstyle": "arc3,rad=0.12",
            },
            bbox={
                "boxstyle": "round,pad=0.28",
                "facecolor": "white",
                "edgecolor": PAL["muted"],
                "linewidth": 0.6,
                "alpha": 0.97,
            },
            zorder=5,
        )

    ax.axvspan(2006, 2008, color=PAL["muted"], alpha=0.12, zorder=0)
    ax.axvline(2006, linestyle="--", color=PAL["muted"], linewidth=0.7, zorder=1)
    ax.axvline(2008, linestyle="--", color=PAL["muted"], linewidth=0.7, zorder=1)
    ymax = max(ys) if ys else 1
    ax.text(
        2007, 9.55, "wwPDB remediation",
        ha="center", va="bottom", fontsize=8, color=PAL["muted"], zorder=5,
    )

    ax.set_xlabel("Deposition year")
    ax.set_ylabel("Mismatches per structure")
    ax.set_xlim(min(xs) - 1.5, max(xs) + 2.2)
    ax.set_ylim(0, max(ymax + 2.6, 10.2))
    _style_axes(ax)

    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
               markersize=7, label=t)
        for t, c in _ERROR_TYPE_COLORS.items()
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=7, frameon=False,
              title="Error class (unlabeled = 1 mismatch)", title_fontsize=7)

    fig.tight_layout()
    out = os.path.join(HERE, "fig2_deposition_year_vs_errors.png")
    _savefig(fig, out)
    plt.close(fig)
    return out


def figure3_signed_volume(vols: Dict[str, np.ndarray]) -> str:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    bins = np.linspace(-4.0, 4.0, 81)
    if vols["d"].size:
        ax.hist(
            vols["d"],
            bins=bins,
            color=PAL["clean"],
            alpha=0.75,
            label=f"D-chirality (n={vols['d'].size:,})",
            edgecolor="none",
        )
    if vols["l"].size:
        ax.hist(
            vols["l"],
            bins=bins,
            color=PAL["error_prone"],
            alpha=0.85,
            label=f"L-chirality (n={vols['l'].size})",
            edgecolor=PAL["navy"],
            linewidth=0.4,
        )
    ax.axvline(0.0, color=PAL["navy"], linewidth=0.9, linestyle="-")
    ax.set_yscale("log")
    ax.set_xlabel(r"Signed C$_\alpha$ tetrahedron volume (Å$^3$)")
    ax.set_ylabel("Residue count")
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    _style_axes(ax)
    fig.tight_layout()
    out = os.path.join(HERE, "fig3_signed_volume_distribution.png")
    _savefig(fig, out)
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

    method_labels: List[str] = []
    method_colors: List[str] = []
    for r in rows:
        m = (r.get("method") or "").upper()
        if "NMR" in m:
            method_labels.append("NMR")
            method_colors.append(PAL["nmr"])
        elif "ELECTRON" in m or "CRYO" in m:
            method_labels.append("Cryo-EM")
            method_colors.append(PAL["cryoem"])
        else:
            method_labels.append("X-ray")
            method_colors.append(PAL["xray"])

    fig, ax = plt.subplots(figsize=(8, 5))
    for method, color in [
        ("X-ray", PAL["xray"]),
        ("NMR", PAL["nmr"]),
        ("Cryo-EM", PAL["cryoem"]),
    ]:
        mask = [lbl == method for lbl in method_labels]
        if not any(mask):
            continue
        ax.scatter(
            mean[mask],
            diff[mask],
            c=color,
            s=55,
            alpha=0.85,
            edgecolor=PAL["navy"],
            linewidth=0.5,
            label=method,
            zorder=3,
        )

    ax.axhline(bias, color=PAL["navy"], linewidth=1.2, linestyle="-", zorder=2)
    ax.axhline(loa_lo, color=PAL["error_prone"], linestyle="--", linewidth=0.9, zorder=2)
    ax.axhline(loa_hi, color=PAL["error_prone"], linestyle="--", linewidth=0.9, zorder=2)
    se = sd / np.sqrt(len(diff)) * 1.96
    ax.fill_between(
        [mean.min() - 0.5, mean.max() + 0.5],
        bias - se,
        bias + se,
        color=PAL["muted"],
        alpha=0.15,
        zorder=1,
    )

    ax.set_xlabel("Mean outlier %")
    ax.set_ylabel("Difference (ChiralFold − wwPDB) %")
    ax.text(
        0.98,
        0.97,
        f"Bias = {bias:.3f}%\nLoA = [{loa_lo:.2f}, {loa_hi:.2f}]%",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color=PAL["navy"],
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": PAL["muted"], "alpha": 0.9},
    )
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    _style_axes(ax)
    fig.tight_layout()
    out = os.path.join(HERE, "fig4_bland_altman.png")
    _savefig(fig, out)
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
    ax_bar.barh(y, rates, color=bar_colors, edgecolor=PAL["navy"], linewidth=0.4, height=0.75)
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(codes, fontsize=9)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Error rate (%)")
    _style_axes(ax_bar)

    vmax = max(int(matrix.max()), 1)
    im = ax_hm.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=vmax)
    ax_hm.set_xticks(range(len(error_types)))
    ax_hm.set_xticklabels(error_types, rotation=35, ha="right", fontsize=9)
    ax_hm.set_yticks(range(len(codes)))
    ax_hm.set_yticklabels(codes, fontsize=9)
    ax_hm.set_xlabel("Error type")

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = int(matrix[i, j])
            if val <= 0:
                continue
            text_color = "white" if val > vmax * 0.55 else PAL["navy"]
            ax_hm.text(
                j,
                i,
                str(val),
                ha="center",
                va="center",
                fontsize=7,
                color=text_color,
                clip_on=True,
            )

    cbar = fig.colorbar(im, ax=ax_hm, fraction=0.035, pad=0.02)
    cbar.set_label("Count", fontsize=9)
    _style_axes(ax_hm)
    fig.tight_layout()
    out = os.path.join(HERE, "fig5_ccd_heatmap.png")
    _savefig(fig, out)
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
    ax.set_xlabel("wwPDB Ramachandran outlier (%)")
    ax.set_ylabel("ChiralFold Ramachandran outlier (%)")
    ax.text(
        0.05,
        0.95,
        f"Spearman ρ = {rho:.2f}\np = {rho_p:.2g}\nn = {len(rows)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color=PAL["navy"],
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": PAL["muted"], "alpha": 0.92},
    )
    ax.legend(loc="lower right", fontsize=9, frameon=False)
    _style_axes(ax)
    fig.tight_layout()
    out = os.path.join(HERE, "fig6_ramachandran_expanded.png")
    _savefig(fig, out)
    plt.close(fig)
    return out


def figure7_af3_resource_benchmark() -> str:
    with open(os.path.join(DATA, "af3_resource_benchmark.json")) as fp:
        data = json.load(fp)

    ops = data["operations"]
    pipeline_ms = ops["full_correct_af3_output_pipeline"]["total_wall_time_s"] * 1000.0
    detection_ms = ops["detection"]["total_wall_time_s"] * 1000.0
    correction_ms = ops["correction"]["total_wall_time_s"] * 1000.0
    residual = data["post_correction_residual_violation_rate_pct"]
    recall = data["aggregate_detection_recall_pct"]
    corrected = data["aggregate_corrected_to_expected_signed_volume_pct"]

    fig, (ax_acc, ax_time) = plt.subplots(1, 2, figsize=(10, 4.2))

    acc_labels = ["Detection recall", "Corrected sign", "Residual violations"]
    acc_vals = [recall, corrected, residual]
    acc_colors = [PAL["clean"], PAL["clean"], PAL["muted"]]
    bars = ax_acc.bar(
        acc_labels,
        acc_vals,
        color=acc_colors,
        edgecolor=PAL["navy"],
        linewidth=0.6,
        width=0.6,
    )
    ax_acc.set_ylim(0, 110)
    ax_acc.set_ylabel("Percent (%)")
    for bar, value in zip(bars, acc_vals):
        ax_acc.text(
            bar.get_x() + bar.get_width() / 2,
            min(value + 3, 104),
            f"{value:.0f}%",
            ha="center",
            va="bottom",
            fontsize=9,
            color=PAL["navy"],
        )
    _style_axes(ax_acc)
    plt.setp(ax_acc.get_xticklabels(), rotation=20, ha="right", fontsize=9)

    time_labels = ["Detection", "Correction", "Full pipeline"]
    time_vals = [detection_ms, correction_ms, pipeline_ms]
    time_colors = [PAL["xray"], PAL["amber"], PAL["navy"]]
    time_bars = ax_time.bar(
        time_labels,
        time_vals,
        color=time_colors,
        edgecolor=PAL["navy"],
        linewidth=0.6,
        width=0.55,
    )
    ax_time.set_ylabel("Wall time (ms)")
    ymax = max(time_vals) * 1.25
    ax_time.set_ylim(0, ymax)
    for bar, value in zip(time_bars, time_vals):
        ax_time.text(
            bar.get_x() + bar.get_width() / 2,
            value + ymax * 0.02,
            f"{value:.1f}",
            ha="center",
            va="bottom",
            fontsize=8,
            color=PAL["navy"],
        )
    _style_axes(ax_time)
    plt.setp(ax_time.get_xticklabels(), rotation=20, ha="right", fontsize=9)

    fig.tight_layout()
    out = os.path.join(HERE, "fig7_af3_resource_benchmark.png")
    _savefig(fig, out)
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
        figure7_af3_resource_benchmark(),
    ]
    for path in outputs:
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
