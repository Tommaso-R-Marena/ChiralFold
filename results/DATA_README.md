# D-Amino Acid Chirality Verification Dataset

## Description

Independent verification of stereochemistry at Cα for all D-amino acid residues in 231 PDB structures. Each residue's chirality was determined by computing the signed volume of the tetrahedron formed by the four Cα substituents (N, C, CB) directly from deposited PDB coordinates.

## Method

```
signed_volume = dot(N − CA, cross(C − CA, CB − CA))
```

- Negative volume → D-chirality (expected for D-labeled residues)
- Positive volume → L-chirality (stereochemistry error)

This verification used **only numpy and raw PDB coordinates**. No ChiralFold code was involved.

## Files

### `d_residue_verification.csv`

1,678 rows, one per D-amino acid residue. Columns:

| Column | Description |
|--------|-------------|
| pdb_id | 4-character PDB accession code |
| chain | Chain identifier |
| resnum | Residue sequence number |
| resname | D-amino acid 3-letter code (DAL, DAR, DAS, etc.) |
| one_letter | One-letter amino acid code |
| has_all_atoms | Whether N, CA, C, CB are all present |
| signed_volume | Signed tetrahedron volume (negative = D, positive = L) |
| chirality | Classification: D, L, flat, or incomplete |
| is_error | True if labeled D but coordinates show L |
| n_xyz | Raw N atom coordinates (x,y,z) |
| ca_xyz | Raw CA atom coordinates (x,y,z) |
| c_xyz | Raw C atom coordinates (x,y,z) |
| cb_xyz | Raw CB atom coordinates (x,y,z) |

### `d_residue_verification_summary.json`

Aggregate statistics including error counts by residue type and PDB structure.

## Results

- **1,677** D-amino acid residues with complete backbone atoms
- **1,656** correctly D (negative signed volume)
- **21** errors (positive signed volume) in **12** PDB structures
- **Error rate: 1.3%**

## Affected Structures

| PDB | Residue | Chain | Position | Signed Volume |
|-----|---------|:-----:|:--------:|:-------------:|
| 1ABI | DPN | I | 56 | +2.49 |
| 1BG0 | DAR | A | 403 | +2.58 |
| 1D7T | DTY | A | 4 | +1.85 |
| 1HHZ | DAL | E | 1 | +2.70 |
| 1KO0 | DLY | A | 542 | +0.12 |
| 1MCB | DHI | P | 3 | +2.60 |
| 1OF6 | DTY | A–H | 1369–1370 | +2.51 to +2.67 |
| 1P52 | DAR | A | 403 | +2.54 |
| 1UHG | DSN | D | 164 | +2.21 |
| 1XT7 | DSG | A | 3 | +2.55 |
| 2AOU | DCY | A | 248 | +2.67 |
| 2ATS | DLY | A | 3001–3003 | +2.56 to +2.59 |

## Reproduction

```bash
python benchmarks/independent_d_residue_verification.py
```

Requires only Python 3.9+ and numpy. PDB files must be in `results/d_survey/`.

## Citation

See `CITATION.cff` in the repository root.

## License

MIT. See repository LICENSE.
