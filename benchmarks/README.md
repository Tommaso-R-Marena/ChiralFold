# ChiralFold Benchmarks

This directory contains the benchmark scripts used to generate the
quantitative results reported in the paper and the README. Every numeric
claim — violation rates, error counts, MolProbity correlations, the
PDB-wide D-residue survey — is reproducible from the scripts below.

## Scripts

| Script | Purpose | Output(s) |
|--------|---------|-----------|
| `independent_d_residue_verification.py` | **Primary independent validation.** Computes the signed Cα tetrahedron volume for every D-amino acid residue in the local PDB cache and flags D-label/L-coordinate mismatches. Uses **only numpy** — no ChiralFold imports. | `results/d_residue_verification.csv`, `results/d_residue_verification_summary.json` |
| `d_residue_survey.py` | Downloads PDB files containing any of the 18 standard D-amino acid CCD codes from RCSB and runs the signed-volume check. Drives the >91% coverage statistic. | PDB files under `results/d_survey/`, intermediate logs |
| `classify_d_errors.py` | Cross-references each mismatch against deposition metadata, COMPND records, and primary literature to assign one of: **Stereochem**, **CCD-Code**, **Mislabel**, **Borderline**. | `results/error_classification.json`, `results/error_table_verified.csv` |
| `childs2025_comparison.py` | Reproduces the 41-sequence ChiralFold vs AlphaFold 3 comparison on the Childs et al. (2025) datasets (DP19/DP9/DP12 + variants). | `results/childs2025_comparison.json` |
| `molprobity_comparison.py` | Benchmarks ChiralFold's auditor against wwPDB/MolProbity validation reports on 31 reference structures (X-ray, NMR, cryo-EM). Reports the canonical Spearman ρ = 0.49, p = 0.006. | `results/molprobity_comparison.json`, `results/molprobity_comparison.png` |
| `bulletproof_verification.py` | Five independent robustness checks confirming the sign convention, the 1KO0 borderline reclassification, the 1OF6 multi-chain consistency, the 1ABI internal control, and a full re-verification of all 16 error structures. | `results/bulletproof_outputs/` |
| `mirror_binder_design.py` | Mirror-image L→D pipeline applied to the p53:MDM2 complex (PDB 1YCR / 3IWY), producing a D-peptide binder candidate retaining the Phe/Trp/Leu hotspot. | `results/3IWY_D_mirror.pdb`, `results/mirror_binder_design.png` |
| `mirror_pipeline_benchmark.py` | L↔D round-trip validation on 5 reference proteins (13,767 atoms total, 0.0 Å coordinate error). | `results/mirror_pipeline_benchmark.png` |
| `folding_quality_benchmark.py` | Conformer-ensemble folding quality (planarity, Ramachandran) on the Childs DP12:MDM2 system. | `results/folding_quality.png` |
| `survey_and_docking_figures.py` | Generates publication-ready figures combining survey statistics and docking scores. | `results/survey_and_docking.png` |
| `pdb_50_set.py` | Downloads the 50-structure reference set used for the MolProbity comparison. | `results/pdb50/` |
| `run_full_benchmark.py` | End-to-end pipeline: runs the chirality benchmark across all 46 sequences (`PURE_D_SEQS` + `DIASTEREOMER_SEQS`) and writes the headline summary. | `results/summary.json`, `results/benchmark_results.png` |
| `test_v22_improvements.py` | Smoke test for the v3.2 planarity-fix and validator improvements. | console output |
| `v3_comprehensive_benchmark.py` | All-in-one v3 benchmark used to verify the full release. | console summary |

## Reproducibility

A complete reproduction of every paper figure and table runs in roughly
**2–6 hours** on a standard laptop. The dominant cost is the PDB downloads
for the D-residue survey (~1,200 files at ~50 KB each). All other steps
are CPU-bound and complete in minutes.

To reproduce from a clean checkout:

```bash
pip install -e ".[dev]"

# 1. Headline ChiralFold vs AF3 benchmark (fast, < 1 min)
python benchmarks/childs2025_comparison.py

# 2. MolProbity comparison (downloads 31 PDBs the first time, ~5 min)
python benchmarks/molprobity_comparison.py

# 3. D-residue survey (~1–2 h dominated by PDB downloads)
python benchmarks/d_residue_survey.py

# 4. Independent verification (numpy only, < 30 s once survey is cached)
python benchmarks/independent_d_residue_verification.py

# 5. Error classification (offline, < 10 s)
python benchmarks/classify_d_errors.py

# 6. Bulletproof verification (offline, < 30 s)
python benchmarks/bulletproof_verification.py

# 7. Full benchmark suite (regenerates summary.json + figures)
python benchmarks/run_full_benchmark.py
```

The independent verification script
(`independent_d_residue_verification.py`) is the primary scientific
validator of the paper's 29-error / 16-structure result. It requires
**only numpy** and the raw PDB coordinate files; no ChiralFold code is
loaded. Reviewers can re-run it directly to confirm the chirality counts.

## Expected runtimes (rough)

| Step | Time | Notes |
|------|------|-------|
| `childs2025_comparison.py` | 30 s | RDKit-only |
| `molprobity_comparison.py` | 3–5 min | first run downloads 31 PDBs |
| `d_residue_survey.py` | 60–120 min | dominated by PDB downloads |
| `independent_d_residue_verification.py` | 10–30 s | local files only |
| `classify_d_errors.py` | 5–10 s | offline metadata join |
| `bulletproof_verification.py` | 10–30 s | offline |
| `run_full_benchmark.py` | 2–4 min | invokes 46-sequence suite |
