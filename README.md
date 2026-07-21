# ChiralFold

**Chirality-correct protein stereochemistry toolkit** — audit any PDB/mmCIF, fix AF3 chirality errors, build D-peptides by construction, and generate exact L↔D mirrors.

![Tests](https://github.com/Tommaso-R-Marena/ChiralFold/actions/workflows/ci.yml/badge.svg?branch=master)
[![PyPI](https://img.shields.io/pypi/v/chiralfold.svg)](https://pypi.org/project/chiralfold/)
[![Hugging Face Space](https://img.shields.io/badge/%F0%9F%A4%97%20Space-Launch%20Web%20UI-0D9488?style=for-the-badge)](https://huggingface.co/spaces/The-Philosopher/ChiralFold-App)
[![Reproduce D-residue errors (5 min)](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Reproduce_PDB_D_Residue_Errors_5min.ipynb)

---

## Start here

| Goal | Link |
|------|------|
| **Use in a browser** | [Hugging Face Space](https://huggingface.co/spaces/The-Philosopher/ChiralFold-App) (in-app **Light / Dark** toggle) |
| **Reproduce D-residue survey errors (&lt;5 min)** | [Colab](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Reproduce_PDB_D_Residue_Errors_5min.ipynb) · `python benchmarks/reproduce_d_residue_errors.py` |
| **Full mmCIF survey** (known errors + live mmCIF-only universe) | [Colab](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Reproduce_mmCIF_D_Residue_Survey.ipynb) |
| **Interactive results dashboard** | [Colab dashboard](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/ChiralFold_Results_Dashboard.ipynb) |
| **Unusual cases + clash-safety** | [Colab](https://colab.research.google.com/github/Tommaso-R-Marena/ChiralFold/blob/master/demos/Demo_Unusual_Cases_Clash_Safety.ipynb) |
| **Web UI (audit / correct / mirror)** | [Hugging Face Space](https://huggingface.co/spaces/The-Philosopher/ChiralFold-App) |

Full reviewer path: [`docs/REVIEWER_COLAB.md`](docs/REVIEWER_COLAB.md).
| **Lean 4 chirality no-go proofs** | [`formal/chirality_nogo/`](formal/chirality_nogo/) |

---

## Install

Works on **Linux, macOS, and Windows** (Python 3.9–3.12).

### pip (recommended)

```bash
pip install -U pip
pip install "chiralfold==3.5.1"
python -c "import chiralfold; print(chiralfold.__version__)"
chiralfold --help
```

### Editable / from Git

```bash
git clone https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install -e .
pip install -e ".[dev]"    # pytest, ruff
pip install -e ".[web]"    # Gradio UI
pip install -e ".[viz]"    # matplotlib/seaborn
```

### conda (RDKit-friendly)

```bash
conda env create -f environment.yml   # or: conda install -c conda-forge rdkit
conda activate chiralfold             # if using environment.yml
pip install "chiralfold==3.5.1"
```

Conda package recipe: [`conda-recipe/`](conda-recipe/).

### Docker

```bash
docker pull ghcr.io/tommaso-r-marena/chiralfold:3.5.1
docker compose up web    # Gradio on http://localhost:7860
# or
docker run --rm -p 7860:7860 ghcr.io/tommaso-r-marena/chiralfold:3.5.1
```

Dockerfiles: `Dockerfile` (CLI/library), `Dockerfile.web` (Gradio).

### Hugging Face Space

Live UI: [The-Philosopher/ChiralFold-App](https://huggingface.co/spaces/The-Philosopher/ChiralFold-App).  
Use the in-app **Light / Dark** toggle (or `?__theme=light` / `?__theme=dark`) — both modes are contrast-checked.  
Source: [`hf_space/`](hf_space/). Local Gradio: `pip install "chiralfold[web]" && chiralfold-web`.

### Colab

Open any notebook under [`demos/`](demos/). Most auto-install `chiralfold==3.5.1` or gemmi/numpy for survey scripts.

Platform notes: [`docs/INSTALL.md`](docs/INSTALL.md) · PyPI publishing: [`docs/PYPI_PUBLISHING.md`](docs/PYPI_PUBLISHING.md).

---

## Reviewer path (&lt;5 minutes, offline)

Reproduce the **29 D-label / L-coordinate mismatches** from the frozen survey CSV — no network beyond clone, numpy only:

```bash
git clone --depth 1 https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install numpy
python benchmarks/reproduce_d_residue_errors.py
```

Expected: `12,573` checkable residues · `29` errors · `16` structures · rate `0.23%`.

---

## What ChiralFold does

```python
from chiralfold import audit_pdb, correct_af3_output, mirror_pdb, ChiralFold

report = audit_pdb("protein.pdb")          # chirality, Rama, clashes, score
correct_af3_output("af3.pdb", "fixed.pdb") # fix inverted stereocenters
mirror_pdb("1YCR.pdb", "1YCR_D.pdb")       # exact L↔D (RMSD 0.0 Å)
pred = ChiralFold().predict("AFWKELDR")     # D-peptide, 0% violations by construction
```

---

## Key results

| Result | Number | Artefact |
|--------|--------|----------|
| PDB-wide D-residue survey (legacy `.pdb`) | **12,573** residues · **29** errors in **16** structures (0.23%) · **4,623** files | `results/d_residue_verification_summary.json` |
| mmCIF known-error re-verify | **29/29** mismatches recovered in **16** structures | `results/mmcif_d_residue_expansion_summary.json` |
| mmCIF-only universe survey | **79** live mmCIF-only D-AA entries · **1** new clear error (**9BC4** DLE, V=+2.82) | `results/mmcif_only_universe_ids.json` |
| Experimental validation | **14/14** non-borderline pass (2 borderline) | `results/experimental_validation_report.json` |
| Ramachandran vs wwPDB (paper) | **n=362** · Spearman **ρ=0.52** · Pearson **r=0.853** | `results/ramachandran_279struct_chainfix_summary.json` |
| AF3 synthetic correction | **100%** detection · **0%** residual · ~37 ms | `results/af3_resource_benchmark.json` |
| Mirror clashscore | **Unchanged** (isometry) | `tests/test_clash_preservation.py` |
| Lean 4 chirality no-go | Distance-only reps cannot recover signed orientation (no `sorry`) | [`formal/chirality_nogo/`](formal/chirality_nogo/) |

MolProbity does **not** flag the D-residue annotation errors (L-only Cα check).

### Why “~245” became a full live survey of ~79

The historical “~245 mmCIF-only” figure counted failed legacy-`.pdb` downloads against a broad DPN/DAS/DAL CCD enumeration. Colab **can** handle a full mmCIF pass. We now discover the **live** mmCIF-only set by RCSB polymer+ligand CCD search + HEAD filter for missing `.pdb` (**79** entries at freeze; all scanned). That closes the format gap for known errors **and** the checkable mmCIF-only universe (`demos/Reproduce_mmCIF_D_Residue_Survey.ipynb`).

### Lean 4 formal proofs

Machine-checked theorems in [`formal/chirality_nogo/`](formal/chirality_nogo/) (Aristotle / Harmonic): reflection preserves pairwise distances but negates oriented volume; any readout factoring only through the distance matrix cannot determine chirality sign. Build with Lake (see that folder’s README). Methods snippet: `formal/chirality_nogo/paper_methods_snippet.tex`.

---

## Unusual cases

| Case | Example | What happens |
|------|---------|--------------|
| Strained / cyclic macrocycle | **1XT7** · **2RMI** | Signed volume still classifies Cα |
| Non-standard ligands (CCD) | **1OF6** · **1BG0** | Coordinates match L; CCD InChI confirms mislabel |
| Ultra-high resolution | **1HHZ** (0.99 Å) | Labeling failure, not density artifact |
| mmCIF-only (new) | **9BC4** DLE | Universe survey; V=+2.82 Å³ |
| Low-res / non-protein Rama | **5M2K** vancomycin | Excluded from protein Rama benchmark by pre-specified rule |

Mirror L↔D is a global isometry → clashscore identical. AF3 correction reflects the violating Cα across the N–C–Cβ plane; residual chirality violations go to **0%**.

---

## Repository map

```
chiralfold/          # Installable Python package
web/                 # Gradio UI (chiralfold-web)
hf_space/            # Hugging Face Space source
demos/               # Colab notebooks
benchmarks/          # Survey / validation / Rama scripts
results/             # Frozen CSV/JSON artefacts (not hand-edited)
formal/chirality_nogo/  # Lean 4 no-go proofs
tests/               # pytest suite
paper/submission/    # Bioinformatics Overleaf package
Overleaf/            # Zipped submission bundle
docs/                # Install + publishing notes
conda-recipe/        # conda-build recipe
Dockerfile*          # CLI and web images
```

Large coordinate caches (`results/mmcif_cache/`, local `d_survey/` downloads) are gitignored; frozen verification CSVs are the source of truth.

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
