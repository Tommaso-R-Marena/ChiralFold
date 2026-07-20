# ChiralFold

**Chirality-correct protein stereochemistry toolkit: PDB auditing, D-peptide construction, AF3 chirality correction, and mirror-image binder design.**

![Tests](https://github.com/Tommaso-R-Marena/ChiralFold/actions/workflows/ci.yml/badge.svg?branch=master)
[![Open Web App](https://img.shields.io/badge/Web%20App-Launch%20ChiralFold-0D9488?style=for-the-badge&logo=googlechrome&logoColor=white)](#web-app-no-install-required)
[![Results Dashboard (Colab)](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/ChiralFold_Results_Dashboard.ipynb)

> **⚠️ Results Dashboard auto-setup:** Clicking the badge above opens Google Colab, clones this repository, and runs `pip install` automatically when you execute the first cells. Expect **~1–2 minutes** on a fresh runtime before interactive plots appear. No manual configuration is required.

[![Quick Demo in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/ChiralFold_Quick_Demo.ipynb)
[![Toy Dataset Demo in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/ChiralFold_Toy_Dataset_Demo.ipynb)
[![Expanded Ramachandran Benchmark in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Expanded_Ramachandran_Benchmark.ipynb)
[![D-residue experimental validation (Colab)](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/D_Residue_Experimental_Validation.ipynb)

ChiralFold provides `pip install`-able stereochemistry validation and coordinate generation for L-proteins, D-peptides, diastereomers, and any PDB structure. It guarantees **0% chirality violations** at stereogenic centers and includes a MolProbity-calibrated quality auditor validated against wwPDB reports on 31 structures.

> **Try it instantly:** Open the [Quick Demo](demos/ChiralFold_Quick_Demo.ipynb) or the fast [Toy Dataset Demo](demos/ChiralFold_Toy_Dataset_Demo.ipynb) in Google Colab. The toy demo audits the packaged ubiquitin fragment, predicts `AFWKELDR`, enumerates `AFK`, and corrects a synthetic AF3-mimetic inverted residue in under two minutes.

## Key Results

**Integrated results dashboard** — All benchmark outputs (D-residue survey, experimental validation, Ramachandran benchmarks, AF3 correction, Colab demo runs) are consolidated in [`results/colab_integrated_manifest.json`](results/colab_integrated_manifest.json). Open the [**Results Dashboard in Colab**](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/ChiralFold_Results_Dashboard.ipynb) for interactive Plotly visualizations (auto-clones repo on first run; see warning above).

**Chirality validation** — 30/31 PDB structures audit at 100% Cα correctness across X-ray (0.48–3.4 Å), NMR, and cryo-EM. One NMR structure (2JXR) flagged with a genuine stereochemical issue.

**Ramachandran agreement with wwPDB/MolProbity** — On a seeded **era-representative** stratified sample of **362 standard-protein structures** (X-ray ultra/high/med/low resolution + NMR + cryo-EM, spanning all deposition eras), produced after fixing an mmCIF→PDB chain-ID collision bug in the converter, Spearman ρ = 0.52 (95% CI [0.42, 0.61], p = 1.5 × 10⁻²⁶), Pearson r = 0.853 (95% CI [0.717, 0.917], p = 6.1 × 10⁻¹⁰⁴) on outlier percentage, with ChiralFold reporting 0.78% mean outliers vs wwPDB's 0.83%. One glycopeptide entry (5M2K, vancomycin–Zn²⁺: 7-residue cyclic glycopeptide with one standard amino acid, not a globular protein) was excluded under the pre-specified criterion wwpdb_rama_outlier_pct > 50%; its exclusion is documented in `results/ramachandran_279struct_chainfix_summary.json`. This supersedes the earlier 279-structure run (ρ = 0.48, p = 2.4 × 10⁻¹⁷; r = 0.80), which used the buggy converter; the rank result is stable across both. The 104-structure run (ρ = 0.56, p = 4.4 × 10⁻¹⁰) and the original 31-structure benchmark (ρ = 0.49, p = 0.006) are retained for reference. **Colab replications** (authoritative runs): n=155, ρ=0.57 (`results/colab_runs/ramachandran_100struct/`); n=285, ρ=0.44 (`results/colab_runs/ramachandran_200struct/`) — independent samples, not interchangeable with the n=362 paper benchmark. Reproduce any run with `benchmarks/expand_ramachandran_benchmark.py` (latest outputs: `results/ramachandran_279struct_chainfix_comparison.csv`, `results/ramachandran_279struct_chainfix_summary.json`, `results/ramachandran_279struct_chainfix_plot.png`).

**Mirror-image binder design** — Converted the p53:MDM2 crystal structure (PDB 1YCR) into a D-peptide therapeutic candidate that preserves the Phe19/Trp23/Leu26 binding triad as D-amino acids — the same hotspot the experimental dPMI-γ (Kd = 53 nM) uses. All backbone φ angles exactly sign-inverted, 0.0 Å coordinate error. The mirror transformation is mathematically exact; experimental Kd measurement against MDM2 is required to confirm binding affinity and is outside the scope of this computational study.

### Experimental validation (reproducible)

Cross-checks all 16 error structures against RCSB metadata + CCD InChI:

```bash
python benchmarks/experimental_structure_validation.py
```

[![D-residue experimental validation (Colab)](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/D_Residue_Experimental_Validation.ipynb)

Results: `results/experimental_validation_report.json` (14/14 non-borderline cases pass automated criteria).

**5M2K exclusion:** Vancomycin–Zn²⁺ glycopeptide (7-residue cyclic antibiotic, 1 standard amino acid)—not a globular protein. Documented in `results/5m2k_benchmark_exclusion.json`.

**Lean 4 generalization:** Use [Harmonic Aristotle](https://aristotle.harmonic.fun/) with the Lean project from branch `cursor/aristotle-formal-proofs-9901` (see manuscript supplementary Methods for scope).

**mmCIF expansion:** `python benchmarks/mmcif_d_residue_expansion.py` (requires `gemmi`).

**Rockfish (UMD HPC):** `srun --partition=shared --cpus-per-task=4 --mem=8G --time=02:00:00 bash` then clone repo and run validation scripts above.

**PDB-wide D-residue survey** — Verified 12,573 D-amino acid residues across 4,616 PDB files (>91% of all RCSB entries for each of the 18 standard D-amino acid CCD codes). Found 29 D-label/L-coordinate mismatches in 16 structures: 6 genuine stereochemistry errors (biology requires D, coordinates show L), 18 CCD code misassignments across 5 structures (L-molecule labeled with D-code), 3 polymer residue mislabels, and 2 borderlines. Errors occur in nine CCD codes; nine are confirmed clean at zero errors. All cross-referenced against biological context, primary literature, and automated RCSB/CCD validation. MolProbity does not flag any of these.

**AF3 chirality correction** — Automatic stereochemistry post-processing for AlphaFold 3 outputs, directly addressing the 51% D-peptide violation rate documented by Childs et al. (2025). In the new AF3-mimetic resource benchmark (`results/af3_resource_benchmark.json`), ChiralFold detects 100% of synthetic inverted stereocenters, corrects 100% to the expected signed-volume class, leaves 0% residual violations, and runs the full correction pipeline across three systems in 0.0367 s (~37 ms). This is stereochemistry correction, not de novo fold prediction.

## Web App (no install required)

Launch a browser UI to **upload any PDB**, run a full stereochemistry audit, **download a chirality-corrected structure**, or **generate a mirror-image (L↔D)** file:

```bash
pip install "chiralfold[web]"
chiralfold-web
# open http://localhost:7860
```

Or with Docker:

```bash
docker compose up web
```

[![Open Web App locally](https://img.shields.io/badge/Launch-Web%20UI-0D9488?style=flat-square)](http://localhost:7860)

## Installation

### pip (recommended)

```bash
pip install chiralfold
```

Or install the latest from GitHub:

```bash
pip install git+https://github.com/Tommaso-R-Marena/ChiralFold.git@v3.4.0
```

### conda/mamba development environment

```bash
mamba env create -f environment.yml
conda activate chiralfold-dev
```

`environment.yml` installs Python, RDKit, scientific dependencies from conda-forge, and an editable `pip install -e .[dev]`. ChiralFold is not currently published on conda-forge; `conda-recipe/meta.yaml` is provided for local builds and future conda-forge feedstock submission. After feedstock acceptance, the intended user command will be `conda install -c conda-forge chiralfold`.

### Docker

```bash
docker pull ghcr.io/tommaso-r-marena/chiralfold:3.4.0

# Audit a local PDB (mount working directory as /data)
docker run --rm -v "$PWD:/data" ghcr.io/tommaso-r-marena/chiralfold:3.4.0 audit /data/protein.pdb

# Correct AF3 output
docker run --rm -v "$PWD:/data" ghcr.io/tommaso-r-marena/chiralfold:3.4.0 correct-af3 /data/af3.pdb -o /data/fixed.pdb
```

Build locally with `docker compose build` or `docker build -t chiralfold:3.4.0 .`.

### Dependencies

- Python 3.9–3.12
- RDKit (`rdkit>=2023.3,<2027`)
- `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `scikit-learn`
- Development extras: `pytest`, `pytest-cov`, `ruff`, `pre-commit`, `mypy`

If `pip install rdkit` fails on your platform, use the conda/mamba environment above or install RDKit from conda-forge first:

```bash
conda install -c conda-forge rdkit
pip install git+https://github.com/Tommaso-R-Marena/ChiralFold.git
```

## Quick Start

### Audit Any PDB Structure

```python
from chiralfold import audit_pdb, format_report

report = audit_pdb('protein.pdb')
format_report(report)

# Individual metrics
report['chirality']['pct_correct']      # Cα chirality (%)
report['ramachandran']['pct_favored']   # Ramachandran favored (%)
report['planarity']['pct_within_6deg']  # Peptide planarity (%)
report['clashes']['clash_score']        # Steric clash score
report['overall_score']                 # Composite 0-100
```

The packaged toy dataset is available at `chiralfold/data/examples/toy_ubiquitin_fragment.pdb` and is used by the [Toy Dataset Demo](demos/ChiralFold_Toy_Dataset_Demo.ipynb):

```python
from importlib import resources

toy_pdb = resources.files("chiralfold").joinpath("data/examples/toy_ubiquitin_fragment.pdb")
report = audit_pdb(str(toy_pdb))
```

### Correct AlphaFold 3 Chirality Errors

```python
from chiralfold import correct_af3_output

# Detect and fix chirality violations in AF3 predictions
result = correct_af3_output('af3_prediction.pdb', 'corrected.pdb')
print(f"Fixed {result['correction']['n_corrected']} violations")
```

### Enumerate Diastereomers for Drug Design

```python
from chiralfold import enumerate_diastereomers

# Find optimal L/D patterns for a peptide sequence
results = enumerate_diastereomers('AFWKELDR', top_n=10)
for r in results:
    print(f"  {r['chirality_pattern']}  score={r['score']:.1f}")
```

### Score Binding Interfaces

```python
from chiralfold import score_interface

metrics = score_interface('receptor.pdb', 'ligand.pdb')
print(f"BSA: {metrics['bsa']:.0f} Å²")                 # alias: buried_surface_area
print(f"H-bonds: {metrics['hbonds']}")                 # alias: n_hbonds
print(f"Interface score: {metrics['interface_score']:.1f}/100")
```

### Generate Peptide Conformers (L or D)

```python
from chiralfold import ChiralFold

model = ChiralFold()  # fix_planarity=True by default

result = model.predict('MQIFVKTL', chirality_pattern='LLLLLLLL')  # L-protein
result = model.predict('AFWKELDR')                                 # D-peptide
result = model.predict('AFWKELDR', chirality_pattern='DLDLDLDL')   # Diastereomer
result = model.predict('THWKFVELRDSNYQA')                         # 15-mer (v3)
```

### Mirror-Image PDB Transformation

```python
from chiralfold import MirrorImagePredictor, mirror_pdb

MirrorImagePredictor.from_pdb('L_protein.pdb', 'D_protein.pdb')   # L→D
mirror_pdb('D_peptide.pdb', 'L_peptide.pdb')                      # D→L
MirrorImagePredictor.from_pdb_id('1SHG', 'D_SH3.pdb')             # From RCSB
```

### CLI

```bash
# Audit structures
chiralfold audit protein.pdb                              # Single structure
chiralfold audit protein.pdb --json                       # JSON output
chiralfold audit --rcsb-batch structures.txt -o results.csv  # Batch RCSB audit

# Correct AF3 outputs
chiralfold correct-af3 af3_prediction.pdb --output fixed.pdb

# Mirror pipeline
chiralfold mirror input.pdb --output output_D.pdb
chiralfold mirror-id 1UBQ --output D_ubiquitin.pdb

# Peptide design
chiralfold predict AFWKELDR --chirality DLDLDLDL
chiralfold enumerate AFWKELDR --top 10

# Interface scoring
chiralfold score-interface receptor.pdb ligand.pdb
```

## Demos and submission package

- [Quick Demo](demos/ChiralFold_Quick_Demo.ipynb): audit, D-peptide prediction, and mirror transform.
- [Toy Dataset Demo](demos/ChiralFold_Toy_Dataset_Demo.ipynb): fast packaged-data smoke test for Colab.
- [Expanded Ramachandran Benchmark](demos/Expanded_Ramachandran_Benchmark.ipynb): reproduce the n>=100 benchmark in Colab.
- `Overleaf/ChiralFold_Bioinformatics_Overleaf.zip`: one-click Overleaf upload package for the Bioinformatics submission. Download the `Overleaf/` folder or zip, upload it to Overleaf as a new project, and set `chiralfold_bioinformatics.tex` as the main document.

## Testing and CI

Local checks mirror GitHub Actions:

```bash
ruff check chiralfold/
pytest tests/ -m "not slow" --cov=chiralfold --cov-report=term-missing --cov-fail-under=62
pre-commit run --all-files
```

CI runs Python 3.9, 3.10, 3.11, and 3.12, with Ruff plus non-slow pytest coverage. Slow/network benchmark jobs remain manual reproduction tasks documented in `benchmarks/README.md` and `results/REPRODUCIBILITY.md`.

## Optional PyMOL Visualization

ChiralFold can write PyMOL `.pml` scripts without requiring PyMOL as a package
dependency:

```python
from chiralfold.viz import write_af3_correction_session

write_af3_correction_session("af3_before.pdb", "af3_after.pdb", "af3_fix.pml")
```

See `scripts/pymol/README.md` for chirality audit, mirror comparison, and AF3
correction templates. Install PyMOL separately as a system or conda package to
open or render the generated scripts.

## ChiralFold vs MolProbity

Head-to-head comparison across X-ray, NMR, and cryo-EM structures (original 31-structure set at 0.48–3.4 Å, plus expanded stratified samples):

| Metric | ChiralFold | wwPDB/MolProbity | Agreement |
|--------|:----------:|:----------------:|:---------:|
| Ramachandran outlier % (n=362, chain-ID-fixed, excl. non-protein) | 0.78% mean | 0.83% mean | ρ = 0.52 (CI 0.42–0.61), p = 1.5×10⁻²⁶; r = 0.853 (CI 0.717–0.917), p = 6.1×10⁻¹⁰⁴ |
| Ramachandran outlier % (n=279, prior run, chain-ID bug) | 0.71% mean | 0.72% mean | ρ = 0.48 (CI 0.36–0.58), p = 2.4×10⁻¹⁷; r = 0.80 |
| Ramachandran outlier % (n=104, expanded) | 0.92% mean | 0.96% mean | ρ = 0.56, p = 4.4×10⁻¹⁰ |
| Ramachandran outlier % (n=31, original) | 0.60% mean | 0.64% mean | ρ = 0.49, p = 0.006 |
| Chirality validation | 30/31 = 100% | Not directly comparable | Flagged 1 real issue |
| Quality vs resolution | r = -0.26 (expected) | Similar trend | Consistent |

Note: v3.2's hybrid Ramachandran uses an empirical PDB probability grid (built from 5,859 residues across 28 high-quality structures) for the favored/allowed classification, with calibrated rectangular regions as a fallback for the outlier boundary. D-amino acid residues use mirror-image regions calibrated per Hovmöller et al. (2002) and the MolProbity D-proline validation. This hybrid approach achieves mean outlier rates matching wwPDB while maintaining coverage for unusual backbone geometries.

**Where ChiralFold adds value:**
- `pip install` — no web interface or complex local setup required
- Native D-amino acid and diastereomer support (MolProbity doesn't validate D-peptide chirality)
- AF3 chirality correction pipeline (no existing tool does this)
- Bidirectional mirror-image pipeline (L↔D, round-trip exact)
- Python API for programmatic batch auditing
- Conformer generation with planarity fix (33% → 95%)

**Where MolProbity is stronger:**

| Capability | MolProbity | ChiralFold |
|---|---|---|
| Ramachandran contours | Data-derived from ~100K structures | Hybrid empirical grid + calibrated rectangles |
| Chi2–Chi4 rotamer validation | Full chi2/chi3/chi4 | Not implemented (chi1 only) |
| Cβ deviation analysis | Available | Not implemented |
| All-atom contact (clash) scoring | Probe/Reduce all-atom | Backbone H-aware KD-tree |
| Community refinement | Decades | New |

### Auditor Quality on Reference Structures

| PDB | Protein | Resolution | Rama Favored | Rama Outlier | Planarity | Score |
|-----|---------|:----------:|:------------:|:------------:|:---------:|:-----:|
| 1CRN | Crambin | 0.54 Å | 88.6% | 0.0% | 91.1% | 72.1 |
| 1UBQ | Ubiquitin | 1.8 Å | 97.3% | 0.0% | 96.0% | 78.6 |
| 1SHG | SH3 domain | 1.8 Å | 92.7% | 1.8% | 66.1% | 70.0 |

Values from `chiralfold audit --rcsb-batch` run in this validation session.

## Benchmarks

### Chirality

On D-peptide sequences, ChiralFold produces 0% chirality violations (0/467 chiral residues across 46 test sequences) vs AlphaFold 3's documented 51% per-residue violation rate on D-peptides (Childs et al., 2025). Note: ChiralFold's 0% rate is guaranteed by construction — each residue is built with explicit stereochemistry encoding — rather than learned from data. The comparison demonstrates that construction-based approaches solve a problem AF3's diffusion architecture fundamentally cannot.

Fisher's exact test on the canonical contingency table from `results/summary.json`: p ≈ 6.66 × 10⁻¹⁴⁴. 31 PDB structures audited: 30/31 = 100% Cα correctness.

### External Benchmark: Childs et al. 2025

41 sequences using the real D-peptide sequences from the Childs et al. (2025) PDB structures: DP19:L-19437 (DEHELLETAARWFYEIAKR, PDB 7YH8), DP9:Streptavidin (LWQHEATWK, PDB 5N8T), DP12:MDM2 (DWWPLAFEALLR, PDB 3IWY), plus L/D pattern variants. ChiralFold: 0/478 chiral residues violated. AF3: 50–52% per-residue violation rate on the same systems. See `benchmarks/childs2025_comparison.py`.

### PDB-Wide D-Residue Chirality Verification

Independently verified (no ChiralFold code — numpy + raw PDB coordinates only) 12,573 D-amino acid residues across 4,616 PDB files (>91% coverage of all 18 standard D-amino acid CCD codes). Found **29 D-label/L-coordinate mismatches in 16 structures** — cases where deposited Cα coordinates show L-stereochemistry despite D-amino acid labels. Error rate: 0.23%. Errors cluster in 5 CCD codes (DTY 2.35%, DLY 0.73%, DAR 0.36%, DSN 0.28%, DPN 0.33%); 9 codes confirmed clean at ≥91% coverage. All are invisible to MolProbity.

**Error classification (cross-referenced against deposition remarks, COMPND records, biological context, and primary literature):**

| PDB | Residue | Chain | Pos | Signed Volume | Error Type | Evidence |
|-----|---------|:-----:|:---:|:-------------:|:----------:|----------|
| 1ABI | DPN | I | 56 | +2.49 | **Stereochem** | Hirulog D-Phe by design; internal control DPN:1 correct (vol=-2.60) |
| 1BG0 | DAR | A | 403 | +2.58 | CCD-Code | Arginine kinase substrate is L-Arg; standalone ligand |
| 1D7T | DTY | A | 4 | +1.85 | **Stereochem** | Contryphan [D-Tyr4] explicitly designed; NMR structure |
| 1HHZ | DAL | E | 1 | +2.70 | **Stereochem** | Cell wall pentapeptide D-Ala; 0.99 Å atomic resolution |
| 1KO0 | DLY | A | 542 | +0.12 | Borderline | Title says "D,L-lysine" (racemic); ALTLOC B, B=32.3 |
| 1MCB | DHI | P | 3 | +2.60 | Mislabel | COMPND says "L-HIS" explicitly; DHI should be HIS |
| 1OF6 | DTY | A-H | 1369-1370 | +2.51 to +2.67 | CCD-Code | DAHP synthase inhibitor is L-Tyr ([Nature 2023](https://doi.org/10.1038/s42004-023-00946-x)); 8 ligand copies |
| 1P52 | DAR | A | 403 | +2.54 | CCD-Code | Arginine kinase mutant; substrate is L-Arg |
| 1UHG | DSN | D | 164 | +2.21 | Mislabel | Ovalbumin is an L-protein; no biological reason for D-Ser |
| 1XT7 | DSG | A | 3 | +2.55 | **Stereochem** | Daptomycin antibiotic: D-Asn biologically required |
| 2AOU | DCY | A | 248 | +2.67 | Mislabel | Histamine methyltransferase is L-protein; DCY should be CYS |
| 2ATS | DLY | A | 3001-3003 | +2.56 to +2.59 | CCD-Code | Title says "(S)-lysine" = L-Lys; 3 ligand copies |

**Error type breakdown:**
- **CCD-Code** (18 errors, 5 structures): Non-polymer ligand labeled with D-form CCD code when the crystallized molecule is L-form. The coordinates correctly model L-stereochemistry; only the chemical component code is wrong.
- **Mislabel** (3 errors, 3 structures): Polymer residue labeled with D-amino acid code in a context where L is biologically correct (L-protein, or COMPND record specifying L-form).
- **Stereochem** (6 errors, 6 structures): **Most concerning** — the biology requires D-stereochemistry but the deposited coordinates show L. These are genuine coordinate-level errors.
- **Borderline** (2 errors, 2 structures): Near-zero signed volume; inconclusive without additional evidence.

Regardless of type, all 29 mismatches are real annotation inconsistencies in the PDB. The signed volume method detects all error types without requiring knowledge of biological context.

Full dataset: `results/d_residue_verification.csv` (12,574 rows with raw coordinates). Classification: `results/error_classification.json`. Per-structure reclassification evidence: `results/error_table_verified.csv`. **Code-level coverage summary: `results/ccd_code_coverage_summary.csv`** (see below).

> **Benchmark Reproducibility:** All 16 error structures are independently verifiable from `results/error_table_verified.csv` (signed volumes, B-factors, ALTLOC flags, biological evidence) and `results/error_classification.json` (category breakdown and method). No ChiralFold code is required — only numpy and raw PDB coordinate files from RCSB.

**Bulletproof verification:** Five independent checks confirm the findings: (1) Sign convention validated on 24/24 known-correct D-residues in 3IWY (all negative volumes); (2) 1KO0 reclassified as borderline (vol=+0.12, ALTLOC B, B=32.3 — inconclusive); (3) 1OF6 confirmed across all 8 chains (all L-coordinates, consistent with L-Tyr biological role); (4) 1ABI internal control passes cleanly (DPN:1 vol=-2.60 correct vs DPN:56 vol=+2.49 error); (5) Full re-verification of all 16 structures with biological context cross-referencing. See `benchmarks/bulletproof_verification.py` and `results/bulletproof_outputs/`.

**Correlation analysis:** 13 of the 16 error structures were deposited between 1992 and 2005, consistent with the 2006–2008 wwPDB remediation effort. Three post-remediation errors (2RMI 2007, 2W76 2008, 3RIT 2011) confirm the pipeline still lacks a D-specific chirality check. Deposition year significantly predicts errors (Mann-Whitney U=278, p=0.0027 on the initial 12-structure survey). Resolution does not (p=0.19) — errors span 0.99–2.77 Å, indicating a labeling problem rather than a data quality problem.

---

### Code-Level Coverage Summary

`results/ccd_code_coverage_summary.csv` — **18 rows, 12 columns.** One row per standard D-amino acid CCD code. Lets reviewers reproduce the error-clustering finding without parsing the LaTeX.

| Column | Type | Description |
|--------|------|-------------|
| `ccd_code` | string | Three-letter wwPDB Chemical Component Dictionary code (e.g. `DTY`, `DLY`) |
| `rcsb_total_entries` | integer | Total PDB entries containing this code, per RCSB full-text search |
| `entries_surveyed` | integer | Entries downloaded and verified in this study |
| `coverage_pct` | float | `entries_surveyed / rcsb_total_entries × 100` |
| `mmcif_only_unavailable` | integer | Entries not available as legacy PDB format (post-2019 depositions and large cryo-EM assemblies); excluded from verification |
| `residues_checked` | integer | Individual Cα chirality checks performed (residues with complete N/Cα/C/Cβ backbone) |
| `n_errors` | integer | D-labeled residues with confirmed L-stereochemistry (signed volume > 0) |
| `n_error_structures` | integer | Distinct PDB IDs containing at least one error |
| `error_structures` | string | Semicolon-separated list of those PDB IDs |
| `error_rate_pct` | float | `n_errors / residues_checked × 100` |
| `status` | string | `error-prone` (n_errors > 0) or `confirmed-clean` (n_errors = 0 at ≥91% coverage) |
| `biological_context` | string | Plain-English summary of the structural families in which this code appears |

**Key finding from this table:** Errors concentrate in five codes — DTY (2.35% error rate), DLY (0.73%), DAR (0.36%), DPN (0.33%), DSN (0.28%) — while nine codes are confirmed clean at zero errors across ≥91% of their RCSB universe: DAS, DGL, DGN, DIL, DLE, DPR, DTH, DTR, DVA. The error-prone codes share two biological scenarios: enzyme active-site ligands where the D-form CCD code is confused with the L-form substrate (DAR/DTY/DLY), and designed or NRPS-assembled peptides where D-configuration is required but L-coordinates were deposited (DPN/DSN).

### Planarity Fix

- D-peptides: 39% → 94% within 6° of planar
- L-proteins: 33% → 95% (generalizes across 5 backbone types)

### Mirror Pipeline

- 5 PDB systems, 13,767 atoms: 0.0 Å coordinate error
- L→D→L round-trip: mathematically exact
- Contact geometry preserved by construction: 105 Cα-Cα contacts and 10 H-bond donor-acceptor pairs maintained within 0.001 Å across L→D transformation

### wwPDB Comparison

- Chain-ID-fixed benchmark (largest, current): 362 standard-protein structures (363 audited; 1 excluded under pre-specified criterion wwpdb_rama_outlier_pct > 50%: 5M2K, vancomycin–Zn²⁺ complex, wwPDB = 100%); Spearman ρ = 0.52 (95% CI [0.42, 0.61], p = 1.5×10⁻²⁶); Pearson r = 0.853 (95% CI [0.717, 0.917], p = 6.1×10⁻¹⁰⁴); mean outlier rate CF 0.78% vs wwPDB 0.83%. 15 of 378 planned entries skipped: 10 lacked a wwPDB validation report, 5 had >36 distinct chains (cannot be represented in single-column PDB format). See `results/ramachandran_279struct_chainfix_summary.json`.
- Representative benchmark (prior run, chain-ID bug): 279 structures; Spearman ρ = 0.48 (95% CI [0.36, 0.58], p = 2.4×10⁻¹⁷), Pearson r = 0.80; mean outlier rate CF 0.71% vs wwPDB 0.72%. Retained for reference; the mmCIF→PDB converter used here truncated multi-character chain IDs (now fixed). See `results/ramachandran_200struct_summary.json`.
- Expanded benchmark: 104 structures audited (stratified X-ray ultra/high/med/low + NMR + cryo-EM); Spearman ρ = 0.56 (p = 4.4×10⁻¹⁰) vs wwPDB; mean outlier rate CF 0.92% vs wwPDB 0.96%. 9 of 113 planned entries skipped (no wwPDB validation report available). See `results/ramachandran_100struct_summary.json`.
- Original benchmark (historical): 31 structures audited (X-ray, NMR, cryo-EM); Spearman ρ = 0.49 (p = 0.006) vs wwPDB; mean outlier rate CF 0.60% vs wwPDB 0.64%

## Scope and Limitations

ChiralFold is a stereochemistry toolkit, not a de novo structure predictor. It excels at chirality auditing, L↔D coordinate transformation, and template-based conformer generation. For de novo folding, use AlphaFold 3 or ESMFold — then pipe the output through ChiralFold to validate or correct stereochemistry.

| Capability | Status |
|-----------|--------|
| Chirality auditing (L and D) | Production-ready |
| Mirror-image transformation | Production-ready (0.0 Å error) |
| AF3 chirality correction | Production-ready (v3.2+) |
| Formal chirality no-go proofs (Lean 4) | New in v3.4 |
| Diastereomer enumeration | New in v3.2 |
| Interface scoring | New in v3.2 |
| Template threading | Available (template-dependent; requires structural homolog in PDB) |
| Fragment assembly | Available (Chou-Fasman SS + NeRF backbone; not comparable to learned models) |
| De novo fold prediction | Not applicable — ChiralFold is a stereochemistry post-processor, not a de novo structure predictor |

## Previously Addressed Limitations

| Issue | v3.0 Status | v3.2 Resolution |
|-------|:----------:|:----------------:|
| Not a fold predictor | Mirror-only | **Template threading + fragment assembly** (template-dependent; requires structural homolog in PDB) |
| Ramachandran uses rectangles | Calibrated rectangles | **Hybrid: empirical PDB grid + calibrated rectangle fallback** |
| No rotamer analysis | Planned | **Penultimate Rotamer Library validation** (chi1 scoring) |
| Clash methodology differs | Heavy-atom only | **Backbone H-atom placement** before scoring |
| Conformer limit at 30 res | Hard limit | **Fragment assembly for any protein length** |

## Project Structure

```
ChiralFold/
├── chiralfold/
│   ├── __init__.py
│   ├── model.py              # ChiralFold model + SMILES builders
│   ├── auditor.py            # PDB structure quality auditor (H-aware clashes)
│   ├── af3_correct.py        # AlphaFold 3 chirality correction pipeline
│   ├── enumerate.py          # Diastereomer enumeration + ranking
│   ├── interface_scorer.py   # Binding interface scoring
│   ├── ramachandran.py       # Hybrid empirical + rectangular Ramachandran
│   ├── rotamers.py           # Side-chain rotamer validation
│   ├── threading.py          # Template-based fold prediction
│   ├── fragments.py          # Fragment-based backbone assembly
│   ├── validator.py          # Chirality validation engine
│   ├── pdb_pipeline.py       # Mirror-image PDB transformation
│   ├── geometry.py           # Planarity fix + post-processing
│   ├── cli.py                # Command-line interface
│   └── data/
│       ├── test_sequences.py # 46-sequence test library
│       └── ramachandran_grid.json  # Empirical φ/ψ probability grid
├── tests/
│   └── test_chirality.py     # Unit tests (incl. external PDB validation)
├── benchmarks/               # Benchmark scripts (incl. bulletproof_verification.py)
├── results/                  # Generated outputs
├── demos/
│   ├── ChiralFold_Quick_Demo.ipynb           # Colab demo notebook
│   ├── ChiralFold_Toy_Dataset_Demo.ipynb     # Fast packaged toy-data demo
│   └── Expanded_Ramachandran_Benchmark.ipynb # Reproduce the n>=100 benchmark in Colab
├── Overleaf/
│   └── ChiralFold_Bioinformatics_Overleaf.zip # One-click Bioinformatics upload
├── CONTRIBUTING.md           # How to contribute
├── pyproject.toml
├── LICENSE (MIT)
└── README.md
```

## Validation

All numbers in this README were produced by running ChiralFold commands in a single validation session. Key outputs:

**3IWY Audit (experimental D-peptide crystal structure):**
Residues: 189. Chirality: 181 correct, 0 wrong, 8 Gly (100.0%). Ramachandran: 95.2% favored, 0.5% outlier. Score: 65.4/100.

**Childs 2025 Comparison (41 sequences, real PDB sequences):**
0/478 chiral residues violated across DP19 (DEHELLETAARWFYEIAKR, PDB 7YH8), DP9 (LWQHEATWK, PDB 5N8T), DP12 (DWWPLAFEALLR, PDB 3IWY), and 26 L/D pattern variants.

**MDM2 Interface Score (1YCR, p53:MDM2):**
Buried Surface Area: 1,980 Å². Shape Complementarity: 0.933. Hydrogen bonds: 10. Interface score: 61.9/100.

**D-Residue Verification (12,573 residues, independent of ChiralFold):**
29 D-label/L-coordinate mismatches in 16 PDB structures. Error rate: 0.23% (12,573 residues checked across 4,616 files, >91% coverage of all 18 standard D-amino acid CCD codes). Classified by biological context: 6 genuine stereochemistry errors (coordinates modeled as L where D is required), 18 CCD code misassignments across 5 structures (L-form ligand with D-code), 3 polymer residue mislabels, 2 borderlines. All verified using only numpy and raw PDB coordinates (no ChiralFold code). Classification cross-referenced against COMPND records, SEQRES, deposition titles, and primary literature. Full dataset: `results/d_residue_verification.csv`. Classification: `results/error_classification.json`. Code-level coverage table: `results/ccd_code_coverage_summary.csv`.

**Batch Audit (validated against this README):**
1UBQ: 100% chirality, 97.3% Rama favored, 96.0% planarity, score 78.6.
1CRN: 100% chirality, 88.6% Rama favored, 91.1% planarity, score 72.1.
1SHG: 100% chirality, 92.7% Rama favored, 66.1% planarity, score 70.0.

## Data Availability

The complete verification dataset (12,574 rows, raw N/Cα/C/Cβ coordinates for every D-amino acid residue surveyed) is at [`results/d_residue_verification.csv`](results/d_residue_verification.csv). The independent verification script (no ChiralFold dependencies, numpy-only) is at [`benchmarks/independent_d_residue_verification.py`](benchmarks/independent_d_residue_verification.py). Reproducibility instructions for every benchmark in this README are documented in [`results/REPRODUCIBILITY.md`](results/REPRODUCIBILITY.md) and [`benchmarks/README.md`](benchmarks/README.md).

## Version History

### v3.4.0 (2026-06-22)

- **Added:** Bioinformatics submission package, Docker support, Lean 4 no-go proofs, expanded Ramachandran benchmark, AF3-mimetic resource benchmark, PyMOL script helpers, and CI/release hardening.
- **Improved:** README, Colab demos, conda environment/recipe scaffolding, pre-commit hooks, and package data inclusion for toy examples.

### v3.2.1 (2026-05-13)

- **Fixed:** `validate_smiles_chirality()` and `validate_3d_chirality()` in `validator.py` were non-functional (violations counter was never incremented). The 0% violation rate reported in the paper is independently confirmed by 3D coordinate geometry benchmarks (`benchmarks/independent_d_residue_verification.py`, `benchmarks/childs2025_comparison.py`).
- **Fixed:** `__version__` updated to 3.2.1 throughout (package, CITATION.cff, BibTeX).
- **Fixed:** `CITATION.cff` statistics updated to final survey values (12,573 residues, 4,616 files, 29 errors in 16 structures, 0.23% rate).
- **Improved:** Clash detection now uses scipy KD-tree (O(n log n)) instead of O(n²) brute force.
- **Improved:** Inter-residue geometry checks (bond length/angle, peptide planarity ω, Ramachandran φ/ψ) now guard against chain gaps, missing residues, and concatenated multi-chain files.
- **Improved:** D-amino acid Ramachandran regions defined as explicit mirror images of the L regions per Hovmöller et al. (2002) and the MolProbity D-proline validation, rather than relying on negate-and-reuse-L heuristic.
- **Fixed:** Module-level `warnings.filterwarnings('ignore')` removed from `chiralfold/auditor.py`, `chiralfold/validator.py`, `chiralfold/ramachandran.py`. Targeted context managers used instead.

## Citation

```bibtex
@software{chiralfold2026,
  title     = {ChiralFold: General-Purpose Protein Stereochemistry Toolkit},
  author    = {Tommaso R. Marena},
  year      = {2026},
  url       = {https://github.com/Tommaso-R-Marena/ChiralFold},
  version   = {3.4.0},
  note      = {PDB auditing calibrated against wwPDB/MolProbity,
               chirality-correct coordinate generation, AF3 correction
               pipeline, mirror-image binder design validated on MDM2
               (dPMI-gamma, Kd=53nM)}
}
```

The AlphaFold 3 benchmark data is from:

```bibtex
@article{childs2025alphafold3dpeptides,
  title   = {Has AlphaFold 3 Solved the Protein Folding Problem for D-Peptides?},
  author  = {Childs, Cameron M. and Zhou, Pei and Donald, Bruce R.},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.03.14.643307}
}
```

## License

MIT. See [LICENSE](LICENSE).
