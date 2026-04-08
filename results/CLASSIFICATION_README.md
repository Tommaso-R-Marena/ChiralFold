# D-Amino Acid Error Classification Dataset

## Overview

This dataset classifies all 21 D-label/L-coordinate mismatches found in a PDB-wide
survey of 1,677 D-amino acid residues across 231 PDB files. Each mismatch was
cross-referenced against the original PDB deposition metadata, COMPND records,
SEQRES, HET/HETATM records, deposition titles, and primary literature to determine
the root cause.

## Error Types

| Type | Count | Structures | Description |
|------|------:|:----------:|-------------|
| **Stereochem** | 4 | 4 | Biology requires D-stereochemistry; deposited coordinates show L. Most concerning. |
| **CCD-Code** | 13 | 4 | Non-polymer ligand labeled with D-form CCD code; crystallized molecule is L-form. |
| **Mislabel** | 3 | 3 | Polymer residue labeled with D-amino acid code in an L-protein or where the COMPND record explicitly specifies L-form. |
| **Borderline** | 1 | 1 | Inconclusive: near-zero signed volume (+0.12), racemic crystallization conditions. |

## Files

### `error_taxonomy.csv` (12 rows, 21 columns)

One row per error-containing PDB structure. Columns:

| Column | Description |
|--------|-------------|
| PDB_ID | 4-character PDB identifier |
| Residue_Code | 3-letter D-amino acid code as deposited (e.g. DTY, DAR, DPN) |
| Chain | Chain identifier(s) |
| Positions | Residue sequence number(s), semicolon-separated |
| N_Errors | Number of individual error residues in this structure |
| Signed_Volume_Range | Signed tetrahedron volume(s) at Cα in ų |
| Resolution | Crystallographic resolution or "NMR" |
| Deposition_Date | PDB deposition date |
| Title | PDB TITLE record |
| In_SEQRES | Whether the D-residue code appears in SEQRES (polymer vs ligand) |
| ATOM_Lines | Number of ATOM records for this residue code (always 0) |
| HETATM_Lines | Number of HETATM records for this residue code |
| Polymer_Status | "In SEQRES (polymer)" or "Not in SEQRES (non-polymer ligand)" |
| Error_Type | Classification: Stereochem, CCD-Code, Mislabel, or Borderline |
| Biological_Context | Why the D-amino acid label was or was not biologically expected |
| Evidence_COMPND | Evidence from PDB COMPND records |
| Evidence_SEQRES | Evidence from PDB SEQRES records |
| Evidence_Literature | Evidence from published literature |
| Evidence_Internal_Control | Within-structure controls (e.g. correct D-residue in same PDB) |
| Correct_Label | What the residue code should be, with rationale |
| Conclusion | Plain-language classification summary |

### `error_taxonomy_expanded.csv` (21 rows, 13 columns)

One row per individual error residue (expanded from multi-residue structures
like 1OF6 ×8 and 2ATS ×3). Columns:

| Column | Description |
|--------|-------------|
| PDB_ID | 4-character PDB identifier |
| Residue_Code | 3-letter D-amino acid code as deposited |
| Chain | Single chain identifier |
| Position | Residue sequence number |
| Signed_Volume | Signed tetrahedron volume at Cα (ų) |
| B_Factor_CA | B-factor of the Cα atom |
| ALTLOC | Alternate conformer indicator (or "No") |
| Has_OXT | Whether the residue has an OXT atom (free carboxylate = standalone ligand) |
| In_SEQRES | Whether the residue code appears in SEQRES |
| ATOM_Lines | Number of ATOM records for this residue code |
| HETATM_Lines | Number of HETATM records for this residue code |
| Error_Type | Classification: Stereochem, CCD-Code, Mislabel, or Borderline |
| Correct_Label | What the residue code should be |

### `error_classification.json`

Machine-readable summary with total counts per error type and methodology notes.

## Reproducing This Dataset

```bash
# From the repository root:
python benchmarks/classify_d_errors.py
```

The script downloads each PDB file from RCSB (cached in `results/pdb_cache/`),
extracts all metadata, computes signed volumes using only numpy, and writes the
three output files. No ChiralFold code is used.

**Dependencies:** Python 3.8+, numpy

## Method

The signed tetrahedron volume at each Cα is:

```
V = (CB − CA) · ((N − CA) × (C − CA))
```

- V > 0 → L-stereochemistry at Cα
- V < 0 → D-stereochemistry at Cα

For correctly built D-amino acids, V should be negative. A positive V indicates
L-stereochemistry despite the D-amino acid label.

## Classification Criteria

Each mismatch is classified by examining:

1. **SEQRES membership** — Is the residue part of the polymer chain or a standalone ligand?
2. **ATOM vs HETATM** — All D-amino acid residues in these structures are HETATM.
3. **Has OXT** — The presence of an OXT atom (free carboxylate) indicates a standalone amino acid ligand, not a peptide-linked residue.
4. **COMPND record** — Does the compound description specify L- or D-form? (e.g. 1MCB: "L-HIS")
5. **Deposition title** — Does the title reveal the intended stereochemistry? (e.g. 2ATS: "(S)-lysine")
6. **Biological context** — Is the D-amino acid biologically expected in this molecule?
7. **Primary literature** — Published papers describing the structure and its ligands.
8. **Internal controls** — Other D-residues in the same structure that are correctly built.

## Citation

If you use this dataset, please cite the ChiralFold repository:

```
Marena, T.R. (2025). ChiralFold: General-Purpose Protein Stereochemistry Toolkit.
https://github.com/Tommaso-R-Marena/ChiralFold
```
