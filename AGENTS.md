# ChiralFold

ChiralFold is a pure-Python scientific toolkit (library + `chiralfold` CLI) for
protein/peptide stereochemistry: PDB quality auditing, D-peptide construction,
AlphaFold 3 chirality correction, diastereomer enumeration, and mirror-image
(L↔D) transformation. There is **no server or web service** — it is a CLI and an
importable package. See `README.md` and `CONTRIBUTING.md` for full usage.

## Cursor Cloud specific instructions

- Python deps are installed into the **system Python** with
  `pip install --break-system-packages -e ".[dev]"` (Ubuntu marks the interpreter
  as externally managed, so the flag is required; no virtualenv is used). This is
  handled by the startup update script — you do not need to reinstall.
- Editable install: source edits under `chiralfold/` take effect immediately with
  no reinstall. Only re-run the install if dependencies in `pyproject.toml` change.
- `rdkit` needs the system libs `libxrender1` and `libxext6`; these are already
  present in the base image.
- Console scripts (`chiralfold`, `pytest`, `ruff`, `mypy`) install to
  `~/.local/bin`, which is **not on PATH by default**. Either run
  `export PATH="$HOME/.local/bin:$PATH"` first, or invoke via
  `python3 -m pytest` / `python3 -m ruff` (there is no `python -m chiralfold`;
  import the package directly, e.g. `python3 -c "from chiralfold import audit_pdb"`).
- Lint: `ruff check chiralfold/` (CI-gated). Pre-commit config is in
  `.pre-commit-config.yaml`.
- Tests: `pytest tests/ -m "not slow"` is the CI gate (fast, ~25s, no network,
  coverage fail-under 75%). Full suite includes `@pytest.mark.slow` tests that
  download structures from RCSB (`files.rcsb.org`) and need network access.
- Optional conda env: `conda env create -f environment.yml` (RDKit from
  conda-forge). Conda-forge recipe skeleton: `conda-recipe/meta.yaml`.
- AF3 resource benchmark (synthetic mimetics; AF3 outputs are not redistributed):
  `python benchmarks/af3_resource_benchmark.py` → `results/af3_resource_benchmark.json`.
- PyMOL session generators live in `chiralfold.viz` / `scripts/pymol/` (optional;
  system PyMOL may fail on Python 3.12 if the distro package still imports `imp`).
- Bioinformatics Overleaf one-click package:
  `Overleaf/ChiralFold_Bioinformatics_Overleaf.zip` (regenerate with
  `scripts/build_overleaf_package.sh`).
- The `enumerate` command scales as 2^N over sequence length (an 8-mer builds 256
  diastereomers) and can run for many minutes. Use short sequences for quick
  checks. Prefer `chiralfold enumerate SEQ --json`.
