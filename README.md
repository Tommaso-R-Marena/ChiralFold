# ChiralFold

**Chirality-correct peptide 3D coordinate generation for D-amino acids and diastereomers.**

ChiralFold is a cheminformatics pipeline that guarantees **0% chirality violations** at stereogenic centers for any combination of L- and D-amino acids. It addresses a documented failure mode in AlphaFold 3, which produces a **51% per-residue chirality violation rate** on D-peptides — equivalent to random assignment ([Childs, Zhou & Donald, 2025](https://doi.org/10.1101/2025.03.14.643307)).

## What ChiralFold Does (and Doesn't Do)

ChiralFold is **not** a de novo fold predictor competing with AlphaFold 3 on structural accuracy. It is a coordinate generation tool that enforces stereochemistry by construction:

- **De novo mode**: Builds peptide SMILES with explicit `[C@H]`/`[C@@H]` stereocenters, then generates 3D conformers via ETKDG + MMFF94
- **Mirror-image mode**: Transforms known L-peptide structures to D-enantiomers via coordinate reflection (RMSD = 0.0 Å)
- **Mixed L/D mode**: Supports per-residue chirality specification for diastereomeric peptides

AF3's diffusion architecture lacks hard stereochemical constraints, so it treats D-residues as noise. ChiralFold guarantees correct stereochemistry but relies on force-field conformer generation for 3D geometry, which has known limitations (see benchmarks below).

| Metric | AlphaFold 3 | ChiralFold |
|--------|:-----------:|:----------:|
| Chirality violation rate | 51% (random) | **0%** (by construction) |
| Mixed L/D support | 50–52% violation | **0%** (per-residue) |
| Bond length RMSD | — | 0.026 Å |
| Bond angle RMSD | — | 1.5–3.0° |
| Statistical test | — | p < 6.7×10⁻¹⁴⁴ |

## Installation

```bash
pip install git+https://github.com/Tommaso-R-Marena/ChiralFold.git
```

Or from source:

```bash
git clone https://github.com/Tommaso-R-Marena/ChiralFold.git
cd ChiralFold
pip install -e .
```

### Dependencies

Python ≥ 3.9, RDKit ≥ 2023.3, NumPy, SciPy, pandas, matplotlib, seaborn, scikit-learn.

## Quick Start

### De Novo Prediction (with planarity fix)

```python
from chiralfold import ChiralFold

model = ChiralFold()  # fix_planarity=True by default

# Pure D-peptide
result = model.predict('AFWKELDR')
print(result['violation_rate'])  # 0.0

# Mixed L/D diastereomer
result = model.predict('AFWKELDR', chirality_pattern='DLDLDLDL')
```

### Mirror-Image PDB Pipeline (highest quality)

The recommended approach when an L-peptide template exists. Produces
crystallographic-quality D-peptide coordinates by reflecting experimental
structures — bypasses all force-field limitations.

```python
from chiralfold import MirrorImagePredictor, mirror_pdb

# From a local PDB file
result = MirrorImagePredictor.from_pdb('L_protein.pdb', 'D_protein.pdb')

# From RCSB by PDB ID
result = MirrorImagePredictor.from_pdb_id('1SHG', 'D_SH3.pdb')
print(result['n_atoms'])      # Number of atoms transformed
print(result['sequence'])     # One-letter sequence

# Low-level: mirror with chain selection
result = mirror_pdb('complex.pdb', 'complex_D.pdb', chains=['A', 'B'])
```

### Mirror-Image from Coordinates

```python
import numpy as np
from chiralfold import ChiralFold

model = ChiralFold()
l_coords = np.random.randn(100, 3)  # L-peptide coordinates
result = model.predict_from_mirror(l_coords, 'AEAAAKEAAA')
print(result['rmsd_to_ideal_mirror'])  # 0.0
```

### SMILES Builders

```python
from chiralfold import d_peptide_smiles, l_peptide_smiles, mixed_peptide_smiles

d_smi = d_peptide_smiles('AFWK')               # All D
l_smi = l_peptide_smiles('AFWK')               # All L
m_smi = mixed_peptide_smiles('AFWK', 'DLDL')   # Mixed
```

### Chirality Validation

```python
from chiralfold import validate_diastereomer

report = validate_diastereomer('AFWKELDR', 'DLDLDLDL')
print(report['valid'])              # True
print(report['smiles_violations'])  # 0
```

### CLI

```bash
chiralfold predict AFWKELDR
chiralfold predict AFWKELDR --chirality DLDLDLDL
chiralfold validate AFWKELDR --chirality DLDLDLDL
chiralfold benchmark
```

## Benchmarks

### Chirality Benchmark (v2)

Benchmarked against AF3 results from [Childs et al. (2025)](https://doi.org/10.1101/2025.03.14.643307):

- **46 sequences** (31 pure D + 15 mixed L/D diastereomers)
- **0/467 chirality violations** (0.00%)
- **217/217 3D geometry checks passed**
- Fisher's exact p = 6.7×10⁻¹⁴⁴ vs AF3's 51% violation rate
- Cohen's h = −1.591 (very large effect size)

### 3D Geometry Quality (v2.1)

Assessed structural validity of generated conformers beyond chirality:

- **Bond length RMSD**: 0.026 Å from ideal values (PDB typical: ~0.02 Å)
- **Bond angle RMSD**: 1.5–3.0° per type (PDB typical: ~1.5°)
- **Peptide planarity**: 39% of ω angles within 6° of trans — this is a known limitation of ETKDG + MMFF94 for peptide backbone generation
- **Conformer diversity**: 0–14 Å Cα RMSD range across ensemble
- **Radius of gyration**: Consistent with compact/random coil scaling

Ground truth: **PDB 3IWY** — crystal structure of the D-peptide dPMI-gamma (DWWPLAFEALLR) bound to MDM2 at 1.9 Å resolution. This is the same DP12:MDM2 system where AF3 showed 50% chirality violations.

### Limitation Status

| Limitation | v2.1 Status | v2.2 Status |
|-----------|:-----------:|:-----------:|
| Peptide bond planarity | 39% within 6° | **94% within 6°** (fixed) |
| Backbone angle bias | MMFF94-limited | Improved via planarity fix |
| No protein context | De novo only | **Mirror-image PDB pipeline** |
| Not a fold predictor | Inherent | Use mirror-image for real folds |

The **mirror-image PDB pipeline** (v2.2) avoids all force-field limitations by reflecting experimental crystal structures. This is the recommended approach when an L-peptide template exists.

### Reproducing Results

```bash
# Chirality benchmark
python benchmarks/run_full_benchmark.py

# 3D geometry quality benchmark
python benchmarks/folding_quality_benchmark.py

# v2.2 improvement verification
python benchmarks/test_v22_improvements.py
```

Results are saved to `results/`.

## Project Structure

```
ChiralFold/
├── chiralfold/                   # Python package
│   ├── __init__.py
│   ├── model.py                  # ChiralFold model + SMILES builders
│   ├── validator.py              # Chirality validation engine
│   ├── pdb_pipeline.py           # Mirror-image PDB transformation
│   ├── geometry.py               # Planarity fix + backbone post-processing
│   ├── cli.py                    # Command-line interface
│   └── data/
│       └── test_sequences.py     # 46-sequence test library
├── tests/
│   └── test_chirality.py         # 31 unit tests
├── benchmarks/
│   ├── run_full_benchmark.py     # Chirality benchmark (6 phases)
│   ├── folding_quality_benchmark.py  # 3D geometry quality
│   └── test_v22_improvements.py  # Planarity + mirror pipeline tests
├── results/                      # Generated outputs
├── pyproject.toml
├── LICENSE (MIT)
└── README.md
```

## Running Tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## Citation

```bibtex
@software{chiralfold2025,
  title     = {ChiralFold: Chirality-Correct Peptide 3D Coordinate Generation},
  author    = {ChiralFold Contributors},
  year      = {2025},
  url       = {https://github.com/Tommaso-R-Marena/ChiralFold},
  version   = {2.2.0},
  note      = {Guarantees 0\% chirality violations for D-peptides and
               diastereomers where AlphaFold 3 shows 51\% violations.\\{}
               Includes mirror-image PDB pipeline for crystallographic-quality
               D-peptide coordinate generation}
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

MIT. See [LICENSE](LICENSE).
