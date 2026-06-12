#!/usr/bin/env python3
"""
ChiralFold journal-submission statistical extensions.

Reads source data from raw GitHub URLs (with local fallback) and writes
all derived statistics to stats_extended.json. Reproducible with seed 42.

Phases:
  1A Ramachandran agreement (Spearman, Pearson, Bland-Altman, OLS, bootstrap)
  1B Deposition year (Mann-Whitney 12 & 16, Cliff's delta, logistic, permutation)
  1C Resolution analysis (Mann-Whitney + Pearson/Spearman vs |signed_volume|)
  1D Signed-volume bimodality (Hartigan dip test + GMM)
  1E CCD error clustering (Fisher exact + chi-squared)
  1F Full expanded MW on 16-structure dataset
"""
from __future__ import annotations

import io
import json
import os
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture
import statsmodels.api as sm

try:
    import diptest
    HAVE_DIPTEST = True
except Exception:
    HAVE_DIPTEST = False

SEED = 42
RAW_BASE = "https://raw.githubusercontent.com/Tommaso-R-Marena/ChiralFold/master"
LOCAL_BASE = Path(__file__).resolve().parent.parent

FILES = {
    "d_residue_verification.csv": "results/d_residue_verification.csv",
    "ccd_code_coverage_summary.csv": "results/ccd_code_coverage_summary.csv",
    "error_table_verified.csv": "results/error_table_verified.csv",
    "error_classification.json": "results/error_classification.json",
    "molprobity_comparison.json": "results/molprobity_comparison.json",
    "error_correlation_data.json": "results/error_correlation_data.json",
    "benchmark_data.csv": "results/benchmark_data.csv",
}


def fetch(name: str) -> bytes:
    """Get a source file, preferring the local repo, then the raw GitHub URL."""
    rel = FILES[name]
    local = LOCAL_BASE / rel
    if local.exists():
        return local.read_bytes()
    url = f"{RAW_BASE}/{rel}"
    with urllib.request.urlopen(url) as r:
        return r.read()


def load_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(fetch(name)))


def load_json(name: str):
    return json.loads(fetch(name).decode("utf-8"))


# ---------- 1A: Ramachandran agreement ----------
def phase_1A(out: dict) -> None:
    print("\n=== 1A: Ramachandran outlier agreement (n=31) ===")
    mp = load_json("molprobity_comparison.json")
    cf = np.array([row["cf_rama_outlier"] for row in mp], dtype=float)
    ww = np.array([row["wwpdb_rama_outlier"] for row in mp], dtype=float)
    n = len(cf)
    print(f"n = {n}")

    rho, p_rho = stats.spearmanr(cf, ww)
    r, p_r = stats.pearsonr(cf, ww)
    print(f"Spearman rho = {rho:.4f}, p = {p_rho:.4g}")
    print(f"Pearson  r   = {r:.4f}, p = {p_r:.4g}")

    diff = cf - ww
    mean_diff = float(np.mean(diff))
    sd_diff = float(np.std(diff, ddof=1))
    loa_lo = mean_diff - 1.96 * sd_diff
    loa_hi = mean_diff + 1.96 * sd_diff
    se_mean = sd_diff / np.sqrt(n)
    tcrit = stats.t.ppf(0.975, n - 1)
    ci_lo = mean_diff - tcrit * se_mean
    ci_hi = mean_diff + tcrit * se_mean
    print(f"Bland-Altman mean diff = {mean_diff:.4f} (95% CI {ci_lo:.4f} .. {ci_hi:.4f})")
    print(f"             LoA       = [{loa_lo:.4f}, {loa_hi:.4f}]")

    X = sm.add_constant(ww)
    ols = sm.OLS(cf, X).fit()
    slope = float(ols.params[1])
    intercept = float(ols.params[0])
    r2 = float(ols.rsquared)
    p_ols = float(ols.f_pvalue)
    print(f"OLS: cf = {intercept:.4f} + {slope:.4f} * ww;  R^2 = {r2:.4f}, p = {p_ols:.4g}")

    rng = np.random.default_rng(SEED)
    boot_rhos = np.empty(10000)
    idx = np.arange(n)
    for i in range(10000):
        s = rng.choice(idx, size=n, replace=True)
        if np.std(cf[s]) == 0 or np.std(ww[s]) == 0:
            boot_rhos[i] = np.nan
            continue
        boot_rhos[i] = stats.spearmanr(cf[s], ww[s]).statistic
    boot_rhos = boot_rhos[~np.isnan(boot_rhos)]
    boot_lo, boot_hi = np.percentile(boot_rhos, [2.5, 97.5])
    print(f"Bootstrap (n=10000) Spearman rho 95% CI = [{boot_lo:.4f}, {boot_hi:.4f}]")

    out.update({
        "spearman_rho": float(rho),
        "spearman_p": float(p_rho),
        "pearson_r": float(r),
        "pearson_p": float(p_r),
        "bland_altman_mean_diff": mean_diff,
        "bland_altman_sd_diff": sd_diff,
        "bland_altman_loa_lower": float(loa_lo),
        "bland_altman_loa_upper": float(loa_hi),
        "bland_altman_ci_lower": float(ci_lo),
        "bland_altman_ci_upper": float(ci_hi),
        "ols_slope": slope,
        "ols_intercept": intercept,
        "ols_r2": r2,
        "ols_p": p_ols,
        "bootstrap_spearman_ci_lower": float(boot_lo),
        "bootstrap_spearman_ci_upper": float(boot_hi),
        "rama_n_paired": int(n),
        "cf_rama_outlier_mean": float(np.mean(cf)),
        "ww_rama_outlier_mean": float(np.mean(ww)),
    })


# ---------- 1B/1F: Deposition year ----------
def phase_1B_1F(out: dict) -> None:
    print("\n=== 1B/1F: Deposition year vs errors ===")
    eco = load_json("error_correlation_data.json")
    rows = [r for r in eco["results"] if r.get("year") is not None]
    # original 12-error survey
    err12 = [r["year"] for r in rows if r["has_error"]]
    clean = [r["year"] for r in rows if not r["has_error"]]
    print(f"Original error N = {len(err12)}; clean N = {len(clean)}")

    # Expanded 16-error set: add deposition years from error_classification.json
    # 2H9E missing from JSON but documented in error_table_verified.csv and paper as 2006
    ec = load_json("error_classification.json")
    extra_years = {pdb: meta["deposited"]
                   for pdb, meta in ec.get("new_in_targeted_expansion", {}).items()}
    extra_years.setdefault("2H9E", 2006)
    print(f"Expanded error structures (post-2005): {extra_years}")
    err16 = err12 + list(extra_years.values())
    err16_sorted = sorted(err16)
    print(f"Full 16-error years (n={len(err16)}): {err16_sorted}")

    err12_arr = np.asarray(err12, dtype=float)
    err16_arr = np.asarray(err16, dtype=float)
    clean_arr = np.asarray(clean, dtype=float)
    u12, p12 = stats.mannwhitneyu(err12_arr, clean_arr, alternative="two-sided")
    u16, p16 = stats.mannwhitneyu(err16_arr, clean_arr, alternative="two-sided")
    print(f"Mann-Whitney U (12 vs clean) = {u12:.1f}, p = {p12:.4g}")
    print(f"Mann-Whitney U (16 vs clean) = {u16:.1f}, p = {p16:.4g}")

    # Cliff's delta for 16-error vs clean
    def cliffs_delta(a, b):
        a = np.asarray(a)
        b = np.asarray(b)
        gt = np.sum(a[:, None] > b[None, :])
        lt = np.sum(a[:, None] < b[None, :])
        return (gt - lt) / (len(a) * len(b))

    cd = cliffs_delta(err16_arr, clean_arr)
    print(f"Cliff's delta (16 vs clean) = {cd:.4f}")

    # Logistic regression on full 16+clean data (year as predictor of has_error)
    years = np.concatenate([err16_arr, clean_arr])
    yvec = np.array([1] * len(err16) + [0] * len(clean), dtype=float)
    X = sm.add_constant(years)
    logit = sm.Logit(yvec, X).fit(disp=0)
    beta_year = float(logit.params[1])
    se_year = float(logit.bse[1])
    or_year = float(np.exp(beta_year))
    ci_lo = float(np.exp(beta_year - 1.96 * se_year))
    ci_hi = float(np.exp(beta_year + 1.96 * se_year))
    p_logit = float(logit.pvalues[1])
    print(f"Logistic OR(year) = {or_year:.4f}  95% CI [{ci_lo:.4f}, {ci_hi:.4f}]  p = {p_logit:.4g}")

    # Permutation test: difference in mean year between error vs clean
    rng = np.random.default_rng(SEED)
    obs_diff = np.mean(err16_arr) - np.mean(clean_arr)
    pool = np.concatenate([err16_arr, clean_arr])
    n_err = len(err16)
    perm = np.empty(10000)
    for i in range(10000):
        rng.shuffle(pool)
        perm[i] = np.mean(pool[:n_err]) - np.mean(pool[n_err:])
    perm_p = float((np.sum(np.abs(perm) >= abs(obs_diff)) + 1) / (10000 + 1))
    print(f"Permutation p (mean-year diff, n=10000) = {perm_p:.4g}; observed diff = {obs_diff:.3f}")

    out.update({
        "deposition_year_mw_u_12struct": float(u12),
        "deposition_year_mw_p_12struct": float(p12),
        "deposition_year_mw_u_16struct": float(u16),
        "deposition_year_mw_p_16struct": float(p16),
        "cliffs_delta": float(cd),
        "logistic_OR": or_year,
        "logistic_95ci_lower": ci_lo,
        "logistic_95ci_upper": ci_hi,
        "logistic_p": p_logit,
        "permutation_p": perm_p,
        "permutation_obs_diff_mean_year": float(obs_diff),
        "err16_years": err16_sorted,
        "clean_years_n": int(len(clean)),
    })


# ---------- 1C: Resolution analysis ----------
def phase_1C(out: dict) -> None:
    print("\n=== 1C: Resolution analysis ===")
    eco = load_json("error_correlation_data.json")
    rows = [r for r in eco["results"] if r.get("resolution") is not None]
    err_res = np.asarray([r["resolution"] for r in rows if r["has_error"]], dtype=float)
    clean_res = np.asarray([r["resolution"] for r in rows if not r["has_error"]], dtype=float)
    u, p = stats.mannwhitneyu(err_res, clean_res, alternative="two-sided")
    print(f"Resolution MW U = {u:.1f}, p = {p:.4g}; err mean = {np.mean(err_res):.3f}, clean = {np.mean(clean_res):.3f}")

    dv = load_csv("d_residue_verification.csv")
    err_pdbs = set(r["pdb_id"] for r in eco["results"] if r["has_error"])
    res_map = {r["pdb_id"]: r["resolution"] for r in eco["results"] if r.get("resolution") is not None}
    dv2 = dv[dv["pdb_id"].isin(res_map.keys())].copy()
    dv2["resolution"] = dv2["pdb_id"].map(res_map)
    dv2["abs_sv"] = dv2["signed_volume"].abs()
    pr_signed = stats.pearsonr(dv2["resolution"], dv2["signed_volume"])
    sr_signed = stats.spearmanr(dv2["resolution"], dv2["signed_volume"])
    pr_abs = stats.pearsonr(dv2["resolution"], dv2["abs_sv"])
    sr_abs = stats.spearmanr(dv2["resolution"], dv2["abs_sv"])
    print(f"Resolution vs signed_volume: Pearson r={pr_signed.statistic:.4f} p={pr_signed.pvalue:.4g}  Spearman rho={sr_signed.statistic:.4f} p={sr_signed.pvalue:.4g}")
    print(f"Resolution vs |signed_volume|: Pearson r={pr_abs.statistic:.4f} p={pr_abs.pvalue:.4g}  Spearman rho={sr_abs.statistic:.4f} p={sr_abs.pvalue:.4g}")

    out.update({
        "resolution_mw_u": float(u),
        "resolution_mw_p": float(p),
        "resolution_err_mean": float(np.mean(err_res)),
        "resolution_clean_mean": float(np.mean(clean_res)),
        "resolution_signed_pearson_r": float(pr_signed.statistic),
        "resolution_signed_pearson_p": float(pr_signed.pvalue),
        "resolution_signed_spearman_rho": float(sr_signed.statistic),
        "resolution_signed_spearman_p": float(sr_signed.pvalue),
        "resolution_abs_pearson_r": float(pr_abs.statistic),
        "resolution_abs_pearson_p": float(pr_abs.pvalue),
        "resolution_abs_spearman_rho": float(sr_abs.statistic),
        "resolution_abs_spearman_p": float(sr_abs.pvalue),
    })


# ---------- 1D: Signed-volume bimodality ----------
def phase_1D(out: dict) -> None:
    print("\n=== 1D: Signed-volume bimodality ===")
    dv = load_csv("d_residue_verification.csv")
    sv = dv["signed_volume"].dropna().to_numpy()
    print(f"N = {len(sv)}")

    if HAVE_DIPTEST:
        dip, p_dip = diptest.diptest(sv)
        print(f"Hartigan dip test: dip = {dip:.6f}, p = {p_dip:.6g}")
    else:
        rng = np.random.default_rng(SEED)
        dip = _dip_statistic(sv)
        nulls = np.array([_dip_statistic(rng.uniform(sv.min(), sv.max(), size=len(sv))) for _ in range(200)])
        p_dip = float(np.mean(nulls >= dip))
        print(f"Permutation-approx dip: dip = {dip:.6f}, p ~ {p_dip:.4f}")

    gmm = GaussianMixture(n_components=2, random_state=SEED, n_init=10).fit(sv.reshape(-1, 1))
    order = np.argsort(gmm.means_.flatten())
    means = gmm.means_.flatten()[order]
    stds = np.sqrt(gmm.covariances_.flatten())[order]
    weights = gmm.weights_[order]
    print(f"GMM means = {means}")
    print(f"GMM stds  = {stds}")
    print(f"GMM weights = {weights}")

    out.update({
        "dip_statistic": float(dip),
        "dip_p": float(p_dip),
        "gmm_component1_mean": float(means[0]),
        "gmm_component1_std": float(stds[0]),
        "gmm_component1_weight": float(weights[0]),
        "gmm_component2_mean": float(means[1]),
        "gmm_component2_std": float(stds[1]),
        "gmm_component2_weight": float(weights[1]),
        "signed_volume_n": int(len(sv)),
    })


def _dip_statistic(x):
    """Crude unimodality dip statistic (max deviation between ECDF and unimodal fit)."""
    x = np.sort(x)
    n = len(x)
    ecdf = np.arange(1, n + 1) / n
    cdf_lin = (x - x[0]) / (x[-1] - x[0] + 1e-12)
    return float(np.max(np.abs(ecdf - cdf_lin)))


# ---------- 1E: CCD code clustering ----------
def phase_1E(out: dict) -> None:
    print("\n=== 1E: CCD code clustering ===")
    cov = load_csv("ccd_code_coverage_summary.csv")
    cov = cov.copy()
    cov["residues_checked"] = cov["residues_checked"].astype(int)
    cov["n_errors"] = cov["n_errors"].astype(int)
    cov["non_error"] = cov["residues_checked"] - cov["n_errors"]
    error_prone = cov[cov["status"] == "error-prone"]
    clean = cov[cov["status"] != "error-prone"]
    a = int(error_prone["n_errors"].sum())
    b = int(error_prone["non_error"].sum())
    c = int(clean["n_errors"].sum())
    d = int(clean["non_error"].sum())
    print(f"2x2 [errors, non-errors] error-prone: {a}, {b}; clean: {c}, {d}")
    odds, p_fisher = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    print(f"Fisher exact OR = {odds:.4f}, p = {p_fisher:.4g}")

    # Chi-squared across all codes
    table = np.column_stack([cov["n_errors"].to_numpy(), cov["non_error"].to_numpy()])
    chi2, p_chi, dof, _ = stats.chi2_contingency(table)
    print(f"Chi-squared = {chi2:.4f}, dof = {dof}, p = {p_chi:.4g}")

    out.update({
        "fisher_exact_ccd_or": float(odds),
        "fisher_exact_ccd_p": float(p_fisher),
        "chisq_ccd_stat": float(chi2),
        "chisq_ccd_dof": int(dof),
        "chisq_ccd_p": float(p_chi),
        "error_prone_codes": error_prone["ccd_code"].tolist(),
        "clean_codes": clean["ccd_code"].tolist(),
    })


def main():
    np.random.seed(SEED)
    out = {"seed": SEED}
    phase_1A(out)
    phase_1B_1F(out)
    phase_1C(out)
    phase_1D(out)
    phase_1E(out)

    # Save
    out_path = Path(__file__).resolve().parent / "stats_extended.json"
    out_path.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nWrote {out_path}")

    # Mirror to the Overleaf folder if it exists in workspace
    overleaf = Path("/home/user/workspace/ChiralFold_Overleaf")
    if overleaf.exists():
        (overleaf / "stats_extended.json").write_text(json.dumps(out, indent=2, default=float))
        (overleaf / "stats_analysis.py").write_text(Path(__file__).read_text())
        print(f"Mirrored to {overleaf}")


if __name__ == "__main__":
    main()
