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
- Lint: `ruff check chiralfold/`. Note it currently reports pre-existing findings
  and is **not** gated by CI — CI (`.github/workflows/ci.yml`) only runs pytest.
- Tests: `pytest tests/ -m "not slow"` is the CI gate (fast, ~20s, no network).
  The full suite (`pytest tests/`) includes `@pytest.mark.slow` tests that
  download structures from RCSB (`files.rcsb.org`) and need network access.
- The `enumerate` command scales as 2^N over sequence length (an 8-mer builds 256
  diastereomers, each with 3D conformers) and can run for many minutes. Use short
  sequences for quick checks. `chiralfold enumerate SEQ --json` prints structured
  results; the plain (non-`--json`) formatter is comparatively quiet.
