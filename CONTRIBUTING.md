# Contributing to ChiralFold

Contributions are welcome. This guide covers common tasks.

## Setup

```bash
git clone https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install -e ".[dev]"
```

## Running Tests

```bash
# Full test suite
pytest tests/ -v

# Quick run (skip slow tests)
pytest tests/ -v -k "not slow"
```

## Adding a PDB Structure to the Benchmark Set

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

3. Run the auditor to verify it works:
   ```python
   from chiralfold import audit_pdb, format_report
   report = audit_pdb('results/pdb50/xxxx.pdb')
   format_report(report)
   ```

4. If adding a D-amino acid structure, also add it to `benchmarks/d_residue_pdb_ids.json`.

## Adding a New D-Peptide Benchmark Sequence

1. Add the sequence to `chiralfold/data/test_sequences.py`:
   ```python
   DIASTEREOMER_SEQS['YOUR_ID'] = {
       'seq': 'YOURSEQUENCE',
       'chirality': 'DLDLDLDLDL',
       'note': 'Description of the test case',
   }
   ```

2. Run the chirality benchmark to verify:
   ```bash
   python benchmarks/run_full_benchmark.py
   ```

## Code Style

- Python 3.9+
- Type hints where practical
- Docstrings on all public functions
- PDB parsing: use raw line parsing (not BioPython) to minimize dependencies

## Submitting Changes

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature`
3. Run tests: `pytest tests/ -v`
4. Commit with a descriptive message
5. Open a pull request

## Reporting Issues

Open an issue on GitHub with:
- ChiralFold version (`chiralfold --version`)
- Minimal reproducing example
- Expected vs actual behavior
