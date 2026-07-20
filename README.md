# ChiralFold

**Chirality-correct protein stereochemistry toolkit** — audit any PDB, fix AF3 chirality errors, build D-peptides by construction, and generate exact L↔D mirrors.

![Tests](https://github.com/Tommaso-R-Marena/ChiralFold/actions/workflows/ci.yml/badge.svg?branch=master)
[![Hugging Face Space](https://img.shields.io/badge/%F0%9F%A4%97%20Space-Launch%20Web%20UI-0D9488?style=for-the-badge)](https://huggingface.co/spaces/The-Philosopher/ChiralFold-App)
[![Reviewer 5-min](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Reviewer_5min_Reproduce.ipynb)
[![Quick Demo](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/ChiralFold_Quick_Demo.ipynb)

---

## Start here

| Goal | Link |
|------|------|
| **Use in a browser** (upload PDB → correct / mirror → download) | [Hugging Face Space](https://huggingface.co/spaces/The-Philosopher/ChiralFold-App) |
| **Reviewer: reproduce PDB D-residue errors in &lt;5 min** | [Colab notebook](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Reviewer_5min_Reproduce.ipynb) · or `python benchmarks/reproduce_d_residue_errors.py` |
| **Interactive results dashboard** | [Colab dashboard](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/ChiralFold_Results_Dashboard.ipynb) ⚠️ auto-installs (~1–2 min) |
| **Unusual cases + clash-safety demo** | [Colab notebook](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Unusual_Cases_and_Clash_Safety.ipynb) |

---

## Install

Works on **Linux, macOS, and Windows** (Python 3.9–3.12).

```bash
# Recommended (latest from GitHub — always current)
pip install "chiralfold @ git+https://github.com/Tommaso-R-Marena/ChiralFold.git"

# Or clone + editable
git clone https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install -e .
```

Core dependencies: `numpy`, `scipy`, `pandas`, `rdkit`.  
Optional: `pip install "chiralfold[web]"` (Gradio UI) · `pip install "chiralfold[viz]"` (matplotlib/seaborn for benchmarks).

**If `rdkit` fails on your platform:**

```bash
conda install -c conda-forge rdkit
pip install "chiralfold @ git+https://github.com/Tommaso-R-Marena/ChiralFold.git"
```

Full platform notes: [`docs/INSTALL.md`](docs/INSTALL.md).

```bash
# Verify
python -c "import chiralfold; print(chiralfold.__version__)"
chiralfold --help
```

---

## Reviewer path (&lt;5 minutes, offline)

Reproduce the **29 D-label / L-coordinate mismatches** from the frozen survey CSV — no network, no ChiralFold install required beyond numpy:

```bash
git clone --depth 1 https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install numpy   # only dependency for this script
python benchmarks/reproduce_d_residue_errors.py
```

Expected output: `12,573` checkable residues · `29` errors · `16` structures · rate `0.23%`.

Or open the [Reviewer 5-min Colab](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Reviewer_5min_Reproduce.ipynb).

---

## What ChiralFold does

```python
from chiralfold import audit_pdb, correct_af3_output, mirror_pdb, ChiralFold

report = audit_pdb("protein.pdb")          # chirality, Rama, clashes, score
correct_af3_output("af3.pdb", "fixed.pdb") # fix inverted stereocenters
mirror_pdb("1YCR.pdb", "1YCR_D.pdb")       # exact L↔D (RMSD 0.0 Å)
pred = ChiralFold().predict("AFWKELDR")     # D-peptide, 0% violations by construction
```

Web UI (local):

```bash
pip install "chiralfold[web]"
chiralfold-web   # http://localhost:7860
```

---

## Key results (at a glance)

| Result | Number | Where |
|--------|--------|-------|
| PDB-wide D-residue survey | **12,573** residues · **29** errors in **16** structures (0.23%) | `results/d_residue_verification_summary.json` |
| Experimental validation | **14/14** non-borderline pass (2 borderline) | `results/experimental_validation_report.json` |
| Ramachandran vs wwPDB (paper) | **n=362** · Spearman **ρ=0.52** · Pearson **r=0.853** | `results/ramachandran_279struct_chainfix_summary.json` |
| AF3 synthetic correction | **100%** detection · **0%** residual · ~37 ms | `results/af3_resource_benchmark.json` |
| Mirror clashscore | **Unchanged** (isometry — distances preserved) | `tests/test_clash_preservation.py` |

MolProbity does **not** flag the D-residue annotation errors (L-only Cα check).

---

## Unusual cases

ChiralFold is tested on structures that break naive assumptions:

| Case | Example | What happens |
|------|---------|--------------|
| **Strained / cyclic macrocycle** | **1XT7** daptomycin · **2RMI** astressin | Signed volume still classifies Cα; Stereochem errors detected |
| **Non-standard ligands (CCD)** | **1OF6** (8× DTY←L-Tyr) · **1BG0** DAR←L-Arg | Coordinates match L; CCD InChI confirms mislabel |
| **Ultra-high resolution** | **1HHZ** (0.99 Å) DAL error | Not a density/resolution artifact |
| **Low-res / non-protein Rama** | **5M2K** vancomycin glycopeptide | Excluded from protein Rama benchmark by pre-specified rule |

See [`demos/Unusual_Cases_and_Clash_Safety.ipynb`](demos/Unusual_Cases_and_Clash_Safety.ipynb) and `results/5m2k_benchmark_exclusion.json`.

### “By construction” does not invent clashes

- **Mirror L↔D** is a global isometry: all pairwise distances are preserved → **clashscore is identical** before and after (`tests/test_clash_preservation.py`).
- **AF3 chirality correction** reflects only the violating Cα across the N–C–Cβ plane, preserving CA–N / CA–C / CA–Cβ bond lengths. Clashscore may change slightly when a bad stereocenter is fixed; residual chirality violations go to **0%**.

---

## Repository map

```
chiralfold/          # Installable Python package
web/                 # Gradio UI (chiralfold-web)
hf_space/            # Hugging Face Space source
demos/               # Colab notebooks (start with Reviewer_5min_Reproduce)
benchmarks/          # Reproducible survey / validation / Rama scripts
results/             # Frozen CSV/JSON artefacts (do not hand-edit)
tests/               # pytest suite
paper/submission/    # Bioinformatics Overleaf package
docs/                # Install notes, navigation
```

---

## Citation

```bibtex
@software{chiralfold2026,
  author = {Marena, Tommaso R.},
  title  = {ChiralFold: Chirality-correct protein stereochemistry toolkit},
  year   = {2026},
  url    = {https://github.com/Tommaso-R-Marena/ChiralFold}
}
```

AF3 D-peptide context: Childs, Zhou & Donald (2025) bioRxiv [10.1101/2025.03.14.643307](https://doi.org/10.1101/2025.03.14.643307).

---

## License

MIT — see [`LICENSE`](LICENSE).
