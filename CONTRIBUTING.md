# Contributing to ChiralFold

Contributions are welcome. This guide covers common tasks.

## Setup

```bash
git clone https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install -e ".[dev]"
```

## Running tests

```bash
# Fast suite (skip network-dependent and long-running benchmarks)
pytest tests/ -v -m "not slow"

# Full suite (includes the external-PDB integration test)
pytest tests/ -v
```

`tests/test_chirality.py::TestValidatorFix` is the regression test for
the v3.2.1 validator bug fix and must remain green on every PR.

## Code style

- Python ≥ 3.9
- Type hints on all public APIs
- Docstrings on every public function (Args, Returns, Raises, Examples
  where helpful)
- PDB parsing: use raw line parsing rather than BioPython to keep
  dependencies minimal
- Lint with `ruff` before pushing:

```bash
ruff check chiralfold/
```

## Reproducing the benchmarks

See `benchmarks/README.md` and `results/REPRODUCIBILITY.md` for
step-by-step instructions. Quick start:

```bash
python benchmarks/childs2025_comparison.py
python benchmarks/molprobity_comparison.py
python benchmarks/independent_d_residue_verification.py
```

## Adding a PDB structure to the benchmark set

1. Download the PDB file:
   ```bash
   curl -sL "https://files.rcsb.org/download/XXXX.pdb" -o results/pdb50/xxxx.pdb
   ```

2. Add the PDB ID to `benchmarks/pdb50_metadata.json`:
   ```json
   {
     "pdb_id": "XXXX",
     "resolution": 1.8,
     "method": "X-RAY DIFFRACTION",
     "downloaded": true
   }
   ```

3. Run the auditor to confirm it works:
   ```python
   from chiralfold import audit_pdb, format_report
   report = audit_pdb('results/pdb50/xxxx.pdb')
   format_report(report)
   ```

4. If adding a D-amino acid structure, also list it in
   `benchmarks/d_residue_pdb_ids.json`.

## Adding a new D-peptide benchmark sequence

1. Add the sequence to `chiralfold/data/test_sequences.py`:
   ```python
   DIASTEREOMER_SEQS['YOUR_ID'] = {
       'seq': 'YOURSEQUENCE',
       'chirality': 'DLDLDLDLDL',
       'note': 'Description of the test case',
   }
   ```

2. Confirm the chirality benchmark stays at 0 violations:
   ```bash
   python benchmarks/run_full_benchmark.py
   ```

## Data artefacts in `results/`

The CSV, JSON, and PNG files under `results/` are canonical, frozen
outputs of the benchmark scripts. **Do not edit them manually.** If a
benchmark needs updating, change the script in `benchmarks/`, rerun it,
and commit the regenerated artefact.

## Submitting changes

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Run tests: `pytest tests/ -v -m "not slow"`
4. Run linter: `ruff check chiralfold/`
5. Commit with a descriptive message
6. Open a pull request against `master`

## Reporting issues

Open an issue on GitHub with:
- ChiralFold version (`chiralfold --version`)
- Minimal reproducing example
- Expected vs actual behaviour

For private correspondence about the paper or scientific findings,
contact <marena@cua.edu>.
