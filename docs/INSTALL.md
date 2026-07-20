# Cross-platform installation

ChiralFold supports **Python 3.9–3.12** on Linux, macOS, and Windows.

## Recommended install

```bash
pip install -U pip
pip install "chiralfold @ git+https://github.com/Tommaso-R-Marena/ChiralFold.git"
python -c "import chiralfold; print(chiralfold.__version__)"
```

Clone + editable (developers / reviewers):

```bash
git clone https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install -e .
pip install -e ".[dev]"    # pytest, ruff
pip install -e ".[web]"    # Gradio UI
pip install -e ".[viz]"    # matplotlib/seaborn for benchmark plots
```

## RDKit notes by platform

| Platform | Notes |
|----------|-------|
| **Linux** | `pip install rdkit` provides wheels for manylinux |
| **macOS** | `pip install rdkit` (arm64 + x86_64 wheels). If it fails: `conda install -c conda-forge rdkit` then pip-install ChiralFold |
| **Windows** | Prefer 64-bit Python from python.org or conda. If `pip install rdkit` fails, use `conda install -c conda-forge rdkit` first |

## Offline reviewer path (numpy only)

```bash
git clone --depth 1 https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install numpy
python benchmarks/reproduce_d_residue_errors.py
```

## Docker

```bash
docker pull ghcr.io/tommaso-r-marena/chiralfold:latest   # after release publish
docker compose up web   # Gradio on :7860
```

## Troubleshooting

- **`ModuleNotFoundError: rdkit`** → install RDKit via conda-forge, then reinstall ChiralFold.
- **Gradio / `HfFolder` import errors** → `pip install "chiralfold[web]"` pins a compatible `huggingface_hub`.
- **Slow Colab first cell** → expected; notebooks clone + pip install (~1–2 min).
