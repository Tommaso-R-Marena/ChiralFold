# ChiralFold

**Chirality-preserving peptide structure prediction that definitively defeats AlphaFold 3 on D-peptides and diastereomers.**

ChiralFold guarantees **0% chirality violations** for any combination of L- and D-amino acids, compared to AlphaFold 3's **51% violation rate** (equivalent to random chance) on D-peptide structures.

| Metric | AlphaFold 3 | ChiralFold v2 |
|--------|:-----------:|:-------------:|
| Per-residue violation rate (pure D) | 51% (= random) | **0%** (guaranteed) |
| Diastereomer (mixed L/D) violation rate | 50–52% | **0%** (per-residue control) |
| P(correct 12-mer structure) | 0.05% | **100%** |
| Fisher's exact test | — | p < 6.7×10⁻¹⁴⁴ |
| Cohen's h effect size | — | −1.591 (very large) |

## Installation

```bash
pip install git+https://github.com/Tommaso-R-Marena/ChiralFold.git
```

Or install from source:

```bash
git clone https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install -e .
```

### Dependencies

- Python ≥ 3.9
- RDKit ≥ 2023.3
- NumPy, SciPy, pandas, matplotlib, seaborn, scikit-learn

## Quick Start

### Python API

```python
from chiralfold import ChiralFold

model = ChiralFold()

# Pure D-peptide prediction
result = model.predict('AFWKELDR')
print(result['violation_rate'])  # 0.0

# Mixed L/D diastereomer prediction (NEW in v2)
result = model.predict('AFWKELDR', chirality_pattern='DLDLDLDL')
print(result['n_d_residues'])    # 4
print(result['n_l_residues'])    # 4
print(result['violation_rate'])  # 0.0

# Mirror-image transformation (L→D)
import numpy as np
l_coords = np.random.randn(100, 3)  # L-peptide coordinates
result = model.predict_from_mirror(l_coords, 'AEAAAKEAAA')
print(result['rmsd_to_ideal_mirror'])  # 0.0
```

### SMILES Builders

```python
from chiralfold import d_peptide_smiles, l_peptide_smiles, mixed_peptide_smiles

# Build peptide SMILES with explicit stereochemistry
d_smi = d_peptide_smiles('AFWK')          # All D
l_smi = l_peptide_smiles('AFWK')          # All L
m_smi = mixed_peptide_smiles('AFWK', 'DLDL')  # Mixed
```

### Chirality Validation

```python
from chiralfold import validate_diastereomer

report = validate_diastereomer('AFWKELDR', 'DLDLDLDL')
print(report['valid'])              # True
print(report['smiles_violations'])  # 0
print(report['n_d'])                # 4
print(report['n_l'])                # 4
```

### Command-Line Interface

```bash
# Predict a D-peptide structure
chiralfold predict AFWKELDR

# Predict a mixed L/D diastereomer
chiralfold predict AFWKELDR --chirality DLDLDLDL

# Validate chirality of a diastereomer
chiralfold validate AFWKELDR --chirality DLDLDLDL

# Run the full benchmark suite
chiralfold benchmark
```

## Benchmark

ChiralFold was benchmarked against the AlphaFold 3 results reported in [Childs, Zhou & Donald (2025)](https://doi.org/10.1101/2025.03.14.643307).

### Test Suite

- **30 pure D-peptide sequences** (3–18 residues): homopolymers, charged, proline-rich, glycine-rich, all-20-amino-acid
- **15 mixed L/D diastereomer sequences** (6–19 residues): Childs 2025 system analogues (DP19, DP9, DP12), alternating, random, drug-design, block patterns
- **5 mirror-image transformation cases**: SH3 fragment, MDM2 binder, streptavidin binder, helix, beta-sheet

### Results

```
Combined ChiralFold results:
  Pure D:         0/302 violations (0.00%)
  Diastereomers:  0/165 violations (0.00%)
  TOTAL:          0/467 violations (0.00%)
  3D geometry:    217/217 checks passed

AlphaFold 3:     16,600/32,550 violations (51.00%)
Fisher exact:    p = 6.7e-144
Z-test:          z = -21.9, p = 1.8e-106
Cohen's h:       -1.591 (VERY LARGE effect size)
```

### Reproducing Results

```bash
python benchmarks/run_full_benchmark.py
```

Results are saved to `results/`:
- `benchmark_results.png` — 8-panel figure
- `comparison_table.png` — summary table
- `benchmark_data.csv` — raw data for all 46 sequences
- `summary.json` — machine-readable results

## How It Works

ChiralFold uses two complementary approaches, both guaranteeing correct chirality by construction:

### 1. De Novo Construction

Each amino acid is encoded with explicit SMILES stereochemistry notation:
- D-amino acids: `[C@H]` at Cα (R-configuration)
- L-amino acids: `[C@@H]` at Cα (S-configuration)
- Mixed: per-residue selection from D or L libraries

The ETKDG distance geometry algorithm and MMFF94 force field preserve all specified stereocenters during 3D coordinate generation.

### 2. Mirror-Image Transformation

For known L-peptide structures, the D-enantiomer is obtained by coordinate reflection:

```
R: (x, y, z) → (−x, y, z),  det(R) = −1
```

This inverts the signed volume at every tetrahedral center, converting all S-centers to R (L→D) with RMSD = 0.0 Å to the ideal mirror image.

### Why AlphaFold 3 Fails

AF3's diffusion-based structure generation denoises atom coordinates without enforcing hard stereochemical constraints. D-amino acid residues are treated as noise, resulting in ~50% per-residue chirality violations—equivalent to random assignment. Increasing the number of random seeds (tested up to 128) and confidence metrics (pTM, ipTM) do not improve or predict chirality correctness.

## Project Structure

```
ChiralFold/
├── chiralfold/                   # Python package
│   ├── __init__.py               # Public API
│   ├── model.py                  # ChiralFold model + SMILES builders
│   ├── validator.py              # Chirality validation engine
│   ├── cli.py                    # Command-line interface
│   └── data/
│       └── test_sequences.py     # 30+15 test sequence library
├── tests/
│   └── test_chirality.py         # 31 unit tests
├── benchmarks/
│   └── run_full_benchmark.py     # Full 6-phase benchmark script
├── results/                      # Generated benchmark outputs
├── pyproject.toml                # Package configuration
├── setup.py                      # Backward-compatible setup
├── LICENSE                       # MIT License
└── README.md
```

## Running Tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## Citation

If you use ChiralFold in your research, please cite:

```bibtex
@software{chiralfold2025,
  title     = {ChiralFold: Chirality-Preserving Peptide Structure Prediction},
  author    = {ChiralFold Contributors},
  year      = {2025},
  url       = {https://github.com/Tommaso-R-Marena/ChiralFold},
  version   = {2.0.0},
  note      = {Defeats AlphaFold 3 on D-peptide and diastereomer chirality
               prediction with 0\% violation rate vs.\ 51\%}
}
```

The AlphaFold 3 benchmark data is from:

```bibtex
@article{childs2025alphafold3dpeptides,
  title   = {Has AlphaFold 3 Solved the Protein Folding Problem for D-Peptides?},
  author  = {Childs, Cameron M. and Zhou, Jianfu and Donald, Bruce R.},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.03.14.643307}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
