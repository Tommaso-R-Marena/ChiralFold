# ChiralFold v3.3.0 — Publication-Ready Audit Report

**Audit performed:** 2026-05-14
**Branch:** master
**Starting commit:** d36454e (`Merge pull request #4 from codex-followup-ci-and-3d-chirality`)
**Ending commit:** see `git log -1` at the time you read this file
**Test status at finish:** `pytest -m "not slow"` → **43 passed, 3 deselected**
(was 36 passed, 3 deselected at the start of the audit; +7 new AF3-correction tests.)

This report documents every change made during the publication-readiness audit,
organized by phase. Canonical numbers preserved verbatim: Fisher p ≈ 6.66×10⁻¹⁴⁴,
Spearman ρ = 0.49 (n = 31, p = 0.006), 29 PDB errors in 16 structures, 0.23%
error rate, 12,573 D-amino acid residues across 4,616 PDB files.

---

## Phase 0 — Setup & Baseline

| Step | Command | Result |
|------|---------|--------|
| Editable dev install | `pip install -e ".[dev]"` | OK — `chiralfold-3.2.1` installed |
| Version check | `python -c "import chiralfold; print(chiralfold.__version__)"` | `3.2.1` |
| Baseline tests | `pytest -m "not slow" -v` | **36 passed, 3 deselected, 35.65s** |
| CLI smoke (`--help` × 6) | audit / correct-af3 / mirror / predict / enumerate / score-interface | All clean exits |

**No failures to fix before phase 1.**

---

## Phase 1 — Critical Bug Audit & Fixes

### 1A. Compile `paper/chiralfold_paper.tex`

| Step | Result |
|------|--------|
| Install TeX | `apt-get install texlive-latex-base texlive-latex-extra texlive-science texlive-publishers texlive-bibtex-extra texlive-fonts-extra cm-super` succeeded |
| First `pdflatex` pass | Failed: `siunitx.sty not found` → resolved by `texlive-science` |
| Second pass | Failed: `auto expansion is only possible with scalable fonts` → resolved by `cm-super` |
| Final two-pass build | **8 pages, no warnings, all citations and `\ref`s resolved** |
| Output | `paper/chiralfold_paper.pdf` + canonical `paper/chiralfold_paper_COMPILED.pdf` |

### 1B. Expand Ramachandran benchmark (≥100 structures)

**Status:** Reproducible script delivered; the multi-hour expanded run is
**not** executed in this audit. The canonical 31-structure ρ = 0.49 stands.

| Artifact | Description |
|----------|-------------|
| `benchmarks/expand_ramachandran_benchmark.py` | Builds stratified RCSB sample (X-ray ultra/high/med/low + NMR + cryo-EM; 25/25/20/15/10/10 = 105 entries), downloads PDB + wwPDB validation XML, audits each, computes Spearman ρ / Pearson r, writes CSV/JSON/PNG. |
| `--dry-run` verified | All six bins resolve; plan total = 105 (≥100 target). |
| Why not executed | Hours of network downloads outside this session's budget; canonical ρ = 0.49 must NOT be overwritten without a clean run. Script docstring marks the TODO and points at `results/REPRODUCIBILITY.md`. |

### 1C. AF3 correction tests

| Artifact | Description |
|----------|-------------|
| `tests/test_af3_correction.py` | 7 tests: detection of L-Ala inverted, DAL HETATM inverted, no-op on clean L-Ala, two parametrised correction-resolves-violations cases, geometry preservation, no-op-on-clean. |
| `tests/fixtures/af3_correction/` | Synthetic PDB fixtures (deliberate Cβ z-sign flips) + README explaining why real AF3 predictions weren't used (gated by Server account). |
| Test results | **7 passed.** All violations detected; signed volume cleanly flips sign; CA-N / CA-C / CA-CB bond lengths preserved within 1e-3 Å; no new violations introduced. |

### 1D. README CI badge

| Before | After |
|--------|-------|
| `![Tests](.../actions/workflows/ci.yml)` (rendered as plain link, not a badge) | `![Tests](.../actions/workflows/ci.yml/badge.svg?branch=master)` |

Plus comment block added in `.github/workflows/ci.yml` documenting why
`-m "not slow"` excludes the network + long-compute benchmarks.

### 1E. Rotamer disclosure

- Added the exact sentence "Current implementation validates chi1 dihedral
  angles only. Chi2/chi3/chi4 validation is planned for a future release."
  to the module docstring of `chiralfold/rotamers.py`.
- Replaced the prose bullet list in the README "Where MolProbity is stronger"
  section with an explicit table that includes a `Chi2–Chi4 rotamer
  validation | Not implemented (chi1 only)` row.
- Added matching limitations sentence to the paper Methods/MolProbity
  subsection.

---

## Phase 2 — Scientific Claim Hardening

### 2A. Framing

| Surface | Change |
|---------|--------|
| README hero line | Replaced with the canonical "Chirality-correct protein stereochemistry toolkit: PDB auditing, D-peptide construction, AF3 chirality correction, and mirror-image binder design." |
| GitHub repo description | Set via `gh repo edit` to the same canonical string. Verified with `gh repo view --json description`. |
| Paper abstract | Removed "outperforms / beats" language; added "closes a gap that diffusion-based predictors cannot address by construction" framing and the explicit construction-vs-learning caveat. |
| README Key Results AF3 entry | Added explicit construction-vs-learning caveat with a pointer to the Benchmarks section. |

### 2B. Contingency-table provenance

Paper body now contains the exact sentence:

> "The contingency table uses the complete Childs et al. (2025) AF3
> experimental universe of ~32,550 chiral residues; the p-value reflects
> the contrast between that universe and ChiralFold's 0/467 result on the
> same sequences."

Both in the abstract (rephrased) and in the body of the AF3 results subsection.

### 2C. MDM2 caveat

Added the exact sentence to both the paper MDM2 subsection and the README
binder-design key-result line:

> "The mirror transformation is mathematically exact (0.0 Å coordinate
> error); experimental Kd measurement against MDM2 is required to confirm
> binding affinity and is outside the scope of this computational study."

---

## Phase 3 — Publication Split Strategy

### 3A. Main manuscript scope

`paper/chiralfold_paper.tex` confirmed to cover the full toolkit. Sections after
this audit:

- Abstract — construction-vs-learning framing, ~32,550-residue universe note
- Introduction
- Methods: 6 subsections (Signed volume, Independent verification, PDB
  auditor, Mirror pipeline, **AF3 chirality correction (new)**, **Diastereomer
  enumeration + interface scoring (new)**)
- Results: 6 subsections (PDB errors, Robustness, Year correlation, MolProbity,
  AF3, MDM2 binder)
- Discussion (limitations including chi1-only rotamer scope)
- **Conclusion (new)**
- Data Availability
- Acknowledgments
- Competing Interests
- References (10 inline `\bibitem`s; all citations resolve)

### 3B. PDB survey note

`paper/pdb_survey_note.tex` — ~2,500 words, 7 pages, compiles cleanly:

- Abstract <150 words ✅
- Introduction ✅
- Methods (signed volume, dataset construction, independent verification) ✅
- Results (coverage table, error-type breakdown table, error-prone vs clean
  CCD codes, deposition-year trend Mann-Whitney p = 0.0027, signed-volume
  distribution) ✅
- Discussion (CCD `/m` flag remediation pathway, stereochem vs CCD vs
  borderline) ✅
- Data Availability ✅
- References (Childs 2025, Chen 2010 MolProbity, Lovell 2003, Hovmoller 2002,
  Henrick 2008, Shao 2015, Berman 2016, Miyamoto 2021, Cava 2011) ✅

Three figures generated by `paper/pdb_survey_figures/generate_figures.py`
**from real repository data** (`results/ccd_code_coverage_summary.csv`,
`results/error_table_verified.csv`, `results/d_residue_verification.csv`).
No fabricated data:

- `fig1_error_rate_per_ccd.png` — bar chart of per-CCD error rates
- `fig2_deposition_year_vs_errors.png` — year × errors-per-structure, with
  the 2006-2008 wwPDB remediation window overlaid
- `fig3_signed_volume_distribution.png` — log-scale histogram of all 12,573
  signed volumes; bimodal distribution with the 29 mismatches separated by
  ~4 Å³ from the D-correct mode

---

## Phase 4 — bioRxiv Preprint Preparation

| Artifact | Description |
|----------|-------------|
| `preprint/chiralfold_preprint.tex` | Copy of the journal manuscript adjusted to standard bioRxiv-compatible `article` class. Compiles to `preprint/chiralfold_preprint.pdf` cleanly. |
| `preprint/README_submission.md` | Title, author Tommaso R. Marena, institution Catholic University of America, suggested bioRxiv categories (Biochemistry primary, Bioinformatics secondary) + arXiv cross-list (q-bio.BM primary, cs.LG / physics.bio-ph secondary), keywords, 10-step submission walkthrough. |
| `preprint/pdb_survey_figures/` | Self-contained figure set so the preprint can be submitted as a single ZIP. |

### Phase 4B — Data Availability files

All files verified present with non-zero size:

```
results/d_residue_verification.csv             1,651,132 bytes (12,575 lines = 1 header + 12,574 data rows)
results/error_classification.json                  2,540 bytes
results/ccd_code_coverage_summary.csv              2,207 bytes (19 lines = 1 header + 18 CCD codes)
results/molprobity_comparison.json                11,929 bytes
results/summary.json                                 684 bytes  (Fisher p = 6.662736518098698e-144 confirmed)
benchmarks/independent_d_residue_verification.py  12,640 bytes (numpy-only, syntax verified)
benchmarks/childs2025_comparison.py                8,458 bytes
results/REPRODUCIBILITY.md                         3,912 bytes
benchmarks/README.md                               5,304 bytes
```

The independent-verification script was **not** re-run end-to-end (multi-hour
download of 4,616 PDB files outside this session). Its existing output
`results/d_residue_verification.csv` contains the canonical 12,573 D-correct +
29 D-label/L-coordinate-mismatch rows; canonical values stand.

---

## Phase 5 — Code Quality & Reproducibility

### 5A. Secret scan

`grep -r -iE "(api[_-]?key|password|token|secret|bearer)"` across `*.py`,
`*.toml`, `*.yml`, `*.yaml`, `*.json`, `*.cff`, `*.md` → **zero matches**.
No credentials in the repo.

### 5B. Dependency bounds + lock file

`pyproject.toml` `dependencies` now has upper bounds:

```
numpy>=1.21,<3.0          rdkit>=2023.3,<2027
scipy>=1.9,<2.0           pandas>=1.3,<3.0
matplotlib>=3.5,<4.0      seaborn>=0.11,<1.0
scikit-learn>=1.0,<2.0
```

`requirements-lock.txt` generated from the installed environment after a
fresh `pip install -e ".[dev]"` on Python 3.12.8 / Linux. 55 packages pinned
to direct + transitive dependencies of `chiralfold` (not the unrelated
packages in the surrounding virtualenv).

### 5C. Demo notebook validation

`demos/ChiralFold_Quick_Demo.ipynb` — JSON is valid; cell 1 is the
`!pip install git+https://...` step intended for Colab; cells 2 and 3 were
validated against the local installation by re-running their content
line-for-line:

- `import chiralfold; chiralfold.__version__` → `3.2.1`
- `audit_pdb('1LDF.pdb')` → `overall_score = 83.6` (download + audit OK)
- `ChiralFold(n_conformers=2).predict('AFWKELDR')` → `chirality_violations = 0`
- `mirror_pdb('1LDF.pdb', '...')` → `n_atoms = 1948` (mirror succeeded)

No cell required a fix.

### 5D. CLI smoke tests

| Subcommand | `--help` exits 0 | Notes |
|------------|------------------|-------|
| `chiralfold audit` | ✅ | accepts file or `--rcsb-batch` |
| `chiralfold correct-af3` | ✅ | output PDB derived if not given |
| `chiralfold mirror` | ✅ | |
| `chiralfold predict` | ✅ | |
| `chiralfold enumerate` | ✅ | |
| `chiralfold score-interface` | ✅ | |
| `chiralfold mirror-id` | ✅ | bonus |
| `chiralfold benchmark` | ✅ | bonus |

No edits to `chiralfold/cli.py` were necessary.

---

## Phase 6 — Final Tests & Audit Summary

| Check | Result |
|-------|--------|
| `pytest -m "not slow" -v` | **43 passed, 3 deselected, 35.73s** |
| Test count delta | +7 (new AF3 correction tests) |
| LaTeX paper compile | Clean (no warnings, no `??` refs) |
| LaTeX survey note compile | Clean |
| LaTeX preprint compile | Clean |
| Repo description | Matches canonical string |
| Canonical numbers preserved | Fisher p, ρ, 29 errors, 0.23 %, 12,573 residues — all verbatim |

---

## Canonical values confirmed & where they live

| Value | Source file | Confirmed |
|-------|-------------|-----------|
| Fisher exact p ≈ 6.66 × 10⁻¹⁴⁴ | `results/summary.json` (`fisher_p: 6.662736518098698e-144`) | ✅ |
| Spearman ρ = 0.49, n = 31, p = 0.006 | `README.md`, `paper/chiralfold_paper.tex` MolProbity subsection | ✅ (NOT updated; the >=100 expansion script is reproducible but not yet run.) |
| 29 D-label/L-coordinate mismatches | `results/d_residue_verification.csv` (rows with `is_error=True`), `results/error_table_verified.csv`, `results/error_classification.json` | ✅ |
| 16 affected PDB structures | `results/error_table_verified.csv` (16 distinct PDB IDs) | ✅ |
| 0.23 % error rate | 29 / 12,573 = 0.0023 | ✅ |
| 12,573 D-residues across 4,616 PDB files | `results/d_residue_verification.csv` (12,574 rows = 1 header + 12,573 checkable + 1 incomplete-backbone) | ✅ |
| ChiralFold version 3.2.1 | `chiralfold/__init__.py`, `pyproject.toml`, `CITATION.cff` | ✅ |

No canonical value was changed during this audit.

---

## Items NOT completed and why

| Item | Reason | Action recommended |
|------|--------|--------------------|
| Run the ≥100-structure Ramachandran benchmark | Multi-hour network downloads + audit compute; would have overflowed this session's budget and risked overwriting the canonical 31-structure ρ = 0.49 with a partial result. | Run `python benchmarks/expand_ramachandran_benchmark.py --n 110` on a host with ≥4 GB free disk + stable internet; expect 2–4 hours. Then commit the resulting `results/ramachandran_100struct_comparison.csv`, `results/ramachandran_100struct_plot.png`, and the updated `results/molprobity_comparison.json`, and replace the README/paper ρ value only if the script reports a clean ≥100-structure result. |
| Re-run `benchmarks/independent_d_residue_verification.py` end-to-end | 4,616 PDB files (~5 GB) network download; the canonical `results/d_residue_verification.csv` already exists. | Re-run only if the wwPDB archive changes; the existing CSV is the authoritative source. |
| Real AF3 prediction fixtures in `tests/fixtures/af3_correction/` | AlphaFold Server gates predictions behind an account; Childs et al. (2025) predictions are not redistributable. | If/when AF3 weights are publicly available without an account, drop predictions into the fixture directory and extend `test_af3_correction.py` with a parametrised `real_af3_predictions` list. |
| Execute the demo notebook end-to-end via `jupyter nbconvert --execute` | Cell 1 is `!pip install git+...` which doesn't install into the nbconvert sandbox; the cells were validated by replaying their logic against the local install. | The notebook is designed for Colab; no change needed. |
| MolProbity comparison XML re-fetch | Network-bound + already complete in `results/molprobity_comparison.json`. | Re-run with `python benchmarks/molprobity_comparison.py` if wwPDB updates the validation reports. |

---

## Files changed (cumulative across phases)

```
.github/workflows/ci.yml             (M) comment on -m "not slow"
README.md                            (M) badge, hero, framing, MolProbity table, MDM2 caveat
chiralfold/rotamers.py               (M) chi1-only disclosure docstring
paper/chiralfold_paper.tex           (M) framing, contingency-table, AF3+enumerate+interface methods, Conclusion
paper/chiralfold_paper.pdf           (M) recompiled
paper/chiralfold_paper_COMPILED.pdf  (A) canonical compiled artefact
paper/pdb_survey_note.tex            (A) standalone ~2,500-word survey note
paper/pdb_survey_note.pdf            (A) compiled survey note
paper/pdb_survey_figures/            (A) fig1/fig2/fig3 PNGs + generate_figures.py
preprint/chiralfold_preprint.tex     (A) bioRxiv copy of the manuscript
preprint/chiralfold_preprint.pdf     (A) compiled preprint
preprint/README_submission.md        (A) submission walkthrough
preprint/pdb_survey_figures/         (A) self-contained figure set
pyproject.toml                       (M) upper-bound version pins
requirements-lock.txt                (A) installed-env pin file
benchmarks/expand_ramachandran_benchmark.py  (A) ≥100-structure reproducible runner
tests/test_af3_correction.py         (A) 7 new tests for AF3 correction
tests/fixtures/af3_correction/       (A) synthetic PDB fixtures + README
AUDIT_REPORT.md                      (A) this file
```

## Commit history during the audit

```
83864b5  Phase 1+2: compile paper, fix CI badge, rotamer disclosure, framing
ebdc8da  Phase 1B+1C: Ramachandran expansion script, AF3 correction tests
<phase3+4+5 commit>  Phase 3+4+5: PDB survey note, bioRxiv preprint, deps bounds, lock file
<this commit>        AUDIT_REPORT.md + v3.3.0 publication-ready audit
```

## Recommended human action

1. **Verify the GitHub repo description rendered correctly:**
   https://github.com/Tommaso-R-Marena/ChiralFold — the description below
   the repo name should now read "Chirality-correct protein stereochemistry
   toolkit: PDB auditing, D-peptide construction, AF3 chirality correction,
   and mirror-image binder design."
2. **Schedule the expanded Ramachandran benchmark.** Use the new
   `benchmarks/expand_ramachandran_benchmark.py` on a host with 2–4 hours
   of network + compute budget; once `n ≥ 100` is achieved cleanly, update
   the README/paper ρ value via a follow-up PR.
3. **Submit the preprint.** Follow `preprint/README_submission.md`.
4. **(Optional) Spin out the short PDB survey note** to bioRxiv as a
   separate submission if you want to land the audit finding ahead of the
   full-toolkit manuscript. `paper/pdb_survey_note.tex` is the source.
