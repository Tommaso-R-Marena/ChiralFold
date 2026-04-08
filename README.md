# ChiralFold

**General-purpose protein stereochemistry toolkit — chirality-correct structure generation, PDB auditing, and mirror-image transformation for any protein.**

ChiralFold started as a D-peptide chirality tool, but v3 generalizes to a complete stereochemistry toolkit that works on regular L-proteins, D-peptides, diastereomers, and any PDB structure. It guarantees **0% chirality violations** at stereogenic centers and provides MolProbity-style quality auditing for any protein.

## What's New in v3

| Feature | v2.x | v3 |
|---------|:----:|:--:|
| D-peptide chirality | 0% violations | 0% violations |
| L-protein support | Limited | **Full** (conformers, audit, mirror) |
| PDB structure auditing | — | **Chirality + geometry + Ramachandran + clashes** |
| Conformer size limit | 10 residues | **30 residues** |
| Mirror pipeline | L→D only | **Bidirectional L↔D** (round-trip exact) |
| Planarity fix | D-peptides | **L and D proteins** (33→95%) |
| Ramachandran scoring | — | **Per-residue + conformer filtering** |
| Clash detection | — | **Full van der Waals overlap check** |

## Installation

```bash
pip install git+https://github.com/Tommaso-R-Marena/ChiralFold.git
```

## Quick Start

### Audit Any PDB Structure

```python
from chiralfold import audit_pdb, format_report

report = audit_pdb('protein.pdb')
format_report(report)

# Access individual metrics
print(report['chirality']['pct_correct'])      # Cα chirality correctness
print(report['ramachandran']['pct_favored'])    # Ramachandran favored %
print(report['planarity']['pct_within_6deg'])   # Peptide planarity
print(report['clashes']['clash_score'])         # Steric clash score
print(report['overall_score'])                  # Composite 0-100
```

### Generate Peptide Conformers (L or D)

```python
from chiralfold import ChiralFold

model = ChiralFold()  # fix_planarity=True by default

# L-protein peptide
result = model.predict('MQIFVKTL', chirality_pattern='LLLLLLLL')

# D-peptide
result = model.predict('AFWKELDR')  # defaults to all-D

# Mixed L/D diastereomer
result = model.predict('AFWKELDR', chirality_pattern='DLDLDLDL')

# Longer peptides now supported (v3: up to 30 residues)
result = model.predict('THWKFVELRDSNYQA')  # 15-mer
```

### Mirror-Image PDB Transformation

```python
from chiralfold import MirrorImagePredictor, mirror_pdb

# L→D transformation
MirrorImagePredictor.from_pdb('L_protein.pdb', 'D_protein.pdb')

# D→L transformation (bidirectional in v3)
mirror_pdb('D_peptide.pdb', 'L_peptide.pdb')

# Download from RCSB and mirror
MirrorImagePredictor.from_pdb_id('1SHG', 'D_SH3.pdb')
```

### CLI

```bash
# Audit any PDB structure
chiralfold audit protein.pdb
chiralfold audit protein.pdb --json

# Mirror a PDB (L↔D)
chiralfold mirror input.pdb --output output_D.pdb
chiralfold mirror-id 1UBQ --output D_ubiquitin.pdb

# Generate peptide structures
chiralfold predict AFWKELDR
chiralfold predict AFWKELDR --chirality LLLLLLLL
chiralfold predict AFWKELDR --chirality DLDLDLDL
```

## PDB Auditor

The auditor validates any protein structure against six quality criteria:

| Check | What It Measures | Ideal |
|-------|-----------------|-------|
| Cα Chirality | Signed volume at each stereocenter | 100% correct |
| Bond Geometry | N-Cα, Cα-C, C-N lengths; N-Cα-C angles vs ideal | RMSD < 0.02 Å |
| Ramachandran | φ/ψ backbone dihedral regions | > 98% favored |
| Peptide Planarity | ω deviation from 180° (trans) | > 95% within 6° |
| Clash Score | Van der Waals overlaps per 1000 atoms | < 10 |
| Overall Score | Weighted composite | 0-100 |

Validated on 5 PDB structures:

| PDB | Protein | Chirality | Ramachandran | Planarity | Score |
|-----|---------|:---------:|:------------:|:---------:|:-----:|
| 1CRN | Crambin (0.54 Å) | 100% | 77% | 91% | 65 |
| 1UBQ | Ubiquitin (1.8 Å) | 100% | 84% | 96% | 70 |
| 1SHG | SH3 domain | 100% | 78% | 66% | 61 |
| 1L2Y | Trp-cage (NMR) | 100% | 72% | 84% | 63 |
| 5HHD | VEGF-A (2.1 Å) | 100% | 81% | 75% | 59 |

## Benchmarks

### Chirality (v2 baseline, still valid)

- 46 sequences, 0/467 chirality violations, p < 6.7×10⁻¹⁴⁴ vs AF3's 51%

### Planarity Fix

- D-peptides: 39% → 94% within 6° of planar
- L-proteins: 33% → 95% within 6° (generalizes to L)

### Mirror Pipeline

- 5 PDB systems, 13,767 atoms: 0.0 Å coordinate error
- L→D→L round-trip: mathematically exact (0.0 Å)

### Extended Conformer Generation

- 12-mer and 15-mer D/L/mixed peptides: 0 chirality violations
- Conformer counts scale adaptively with length

## Limitations

| Issue | Status |
|-------|--------|
| Not a fold predictor | Use mirror-image for real folds |
| MMFF94 backbone bias | Planarity fix addresses the worst issue |
| Conformer limit | Extended to 30 residues (was 10) |
| No protein context | Mirror pipeline inherits crystal packing |
| Clash score high for de novo | Expected for force-field conformers |

## Project Structure

```
ChiralFold/
├── chiralfold/
│   ├── __init__.py
│   ├── model.py              # ChiralFold model + SMILES builders
│   ├── auditor.py            # PDB structure quality auditor (v3)
│   ├── ramachandran.py       # Ramachandran scoring + filtering (v3)
│   ├── validator.py          # Chirality validation engine
│   ├── pdb_pipeline.py       # Mirror-image PDB transformation
│   ├── geometry.py           # Planarity fix + post-processing
│   ├── cli.py                # Command-line interface
│   └── data/
│       └── test_sequences.py # 46-sequence test library
├── tests/
│   └── test_chirality.py     # 31 unit tests
├── benchmarks/               # Benchmark scripts
├── results/                  # Generated outputs
├── pyproject.toml
├── LICENSE (MIT)
└── README.md
```

## Citation

```bibtex
@software{chiralfold2025,
  title     = {ChiralFold: General-Purpose Protein Stereochemistry Toolkit},
  author    = {ChiralFold Contributors},
  year      = {2025},
  url       = {https://github.com/Tommaso-R-Marena/ChiralFold},
  version   = {3.0.0},
  note      = {PDB auditing, chirality-correct coordinate generation,
               and mirror-image transformation for any protein}
}
```

## License

MIT. See [LICENSE](LICENSE).
