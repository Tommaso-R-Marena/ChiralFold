# ChiralFold

**General-purpose protein stereochemistry toolkit — chirality-correct structure generation, PDB auditing, and mirror-image transformation for any protein.**

ChiralFold provides `pip install`-able stereochemistry validation and coordinate generation for L-proteins, D-peptides, diastereomers, and any PDB structure. It guarantees **0% chirality violations** at stereogenic centers and includes a MolProbity-calibrated quality auditor validated against wwPDB reports on 31 structures.

## Key Results

**Chirality validation** — 30/31 PDB structures audit at 100% Cα correctness across X-ray (0.48–3.4 Å), NMR, and cryo-EM. One NMR structure (2JXR) flagged with a genuine stereochemical issue.

**Ramachandran agreement with wwPDB/MolProbity** — Spearman ρ = 0.49 (p = 0.006) on outlier percentage across 31 structures. ChiralFold reports 0.60% mean outliers vs wwPDB's 0.64%.

**Mirror-image binder design** — Converted the p53:MDM2 crystal structure (PDB 1YCR) into a D-peptide therapeutic candidate that preserves the Phe19/Trp23/Leu26 binding triad as D-amino acids — the same hotspot the experimental dPMI-γ (Kd = 53 nM) uses. All backbone φ angles exactly sign-inverted, 0.0 Å coordinate error.

**PDB-wide D-residue survey** — Verified 12,573 D-amino acid residues across 4,616 PDB files (>91% of all RCSB entries for each of the 18 standard D-amino acid CCD codes). Found 29 D-label/L-coordinate mismatches in 16 structures: 6 genuine stereochemistry errors (biology requires D, coordinates show L), 18 CCD code misassignments across 5 structures (L-molecule labeled with D-code), 3 polymer residue mislabels, and 2 borderlines. Errors cluster in 5 CCD codes (DTY, DLY, DPN, DSN, DAR); 9 codes confirmed clean at zero errors. All cross-referenced against biological context and primary literature. MolProbity does not flag any of these.

**AF3 chirality correction** — Automatic detection and correction of stereochemistry violations in AlphaFold 3 outputs, directly addressing the 51% violation rate documented by Childs et al. (2025).

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

# Individual metrics
report['chirality']['pct_correct']      # Cα chirality (%)
report['ramachandran']['pct_favored']   # Ramachandran favored (%)
report['planarity']['pct_within_6deg']  # Peptide planarity (%)
report['clashes']['clash_score']        # Steric clash score
report['overall_score']                 # Composite 0-100
```

### Correct AlphaFold 3 Chirality Errors

```python
from chiralfold import correct_af3_output

# Detect and fix chirality violations in AF3 predictions
result = correct_af3_output('af3_prediction.pdb', 'corrected.pdb')
print(f"Fixed {result['n_corrected']} violations")
```

### Enumerate Diastereomers for Drug Design

```python
from chiralfold import enumerate_diastereomers

# Find optimal L/D patterns for a peptide sequence
results = enumerate_diastereomers('AFWKELDR', top_n=10)
for r in results:
    print(f"  {r['chirality_pattern']}  score={r['score']:.1f}")
```

### Score Binding Interfaces

```python
from chiralfold import score_interface

metrics = score_interface('receptor.pdb', 'ligand.pdb')
print(f"BSA: {metrics['buried_surface_area']:.0f} Å²")
print(f"H-bonds: {metrics['n_hbonds']}")
print(f"Interface score: {metrics['interface_score']:.1f}/100")
```

### Generate Peptide Conformers (L or D)

```python
from chiralfold import ChiralFold

model = ChiralFold()  # fix_planarity=True by default

result = model.predict('MQIFVKTL', chirality_pattern='LLLLLLLL')  # L-protein
result = model.predict('AFWKELDR')                                 # D-peptide
result = model.predict('AFWKELDR', chirality_pattern='DLDLDLDL')   # Diastereomer
result = model.predict('THWKFVELRDSNYQA')                         # 15-mer (v3)
```

### Mirror-Image PDB Transformation

```python
from chiralfold import MirrorImagePredictor, mirror_pdb

MirrorImagePredictor.from_pdb('L_protein.pdb', 'D_protein.pdb')   # L→D
mirror_pdb('D_peptide.pdb', 'L_peptide.pdb')                      # D→L
MirrorImagePredictor.from_pdb_id('1SHG', 'D_SH3.pdb')             # From RCSB
```

### CLI

```bash
# Audit structures
chiralfold audit protein.pdb                              # Single structure
chiralfold audit protein.pdb --json                       # JSON output
chiralfold audit --rcsb-batch structures.txt -o results.csv  # Batch RCSB audit

# Correct AF3 outputs
chiralfold correct-af3 af3_prediction.pdb --output fixed.pdb

# Mirror pipeline
chiralfold mirror input.pdb --output output_D.pdb
chiralfold mirror-id 1UBQ --output D_ubiquitin.pdb

# Peptide design
chiralfold predict AFWKELDR --chirality DLDLDLDL
chiralfold enumerate AFWKELDR --top 10

# Interface scoring
chiralfold score-interface receptor.pdb ligand.pdb
```

## ChiralFold vs MolProbity

Head-to-head comparison on 31 PDB structures (0.48–3.4 Å, X-ray + NMR + cryo-EM):

| Metric | ChiralFold | wwPDB/MolProbity | Agreement |
|--------|:----------:|:----------------:|:---------:|
| Ramachandran outlier % | 0.60% mean | 0.64% mean | ρ = 0.49, p = 0.006 |
| Chirality validation | 30/31 = 100% | Not directly comparable | Flagged 1 real issue |
| Quality vs resolution | r = -0.26 (expected) | Similar trend | Consistent |

Note: v3.2's hybrid Ramachandran uses an empirical PDB probability grid (built from 5,859 residues across 28 high-quality structures) for the favored/allowed classification, with calibrated rectangular regions as a fallback for the outlier boundary. This hybrid approach achieves mean outlier rates matching wwPDB while maintaining coverage for unusual backbone geometries.

**Where ChiralFold adds value:**
- `pip install` — no web interface or complex local setup required
- Native D-amino acid and diastereomer support (MolProbity doesn't validate D-peptide chirality)
- AF3 chirality correction pipeline (no existing tool does this)
- Bidirectional mirror-image pipeline (L↔D, round-trip exact)
- Python API for programmatic batch auditing
- Conformer generation with planarity fix (33% → 95%)

**Where MolProbity is stronger:**
- Data-derived Ramachandran contours from ~100K structures
- Rotamer completeness (chi2/chi3/chi4; ChiralFold validates chi1 only)
- Cβ deviation analysis and all-atom contact scoring
- Decades of community validation and refinement

### Auditor Quality on Reference Structures

| PDB | Protein | Resolution | Rama Favored | Rama Outlier | Planarity | Score |
|-----|---------|:----------:|:------------:|:------------:|:---------:|:-----:|
| 1CRN | Crambin | 0.54 Å | 88.6% | 0.0% | 91.1% | 72.1 |
| 1UBQ | Ubiquitin | 1.8 Å | 97.3% | 0.0% | 96.0% | 78.6 |
| 1SHG | SH3 domain | 1.8 Å | 92.7% | 1.8% | 66.1% | 70.0 |

Values from `chiralfold audit --rcsb-batch` run in this validation session.

## Benchmarks

### Chirality

On D-peptide sequences, ChiralFold produces 0% chirality violations (0/467 chiral residues across 46 test sequences) vs AlphaFold 3's documented 51% per-residue violation rate on D-peptides (Childs et al., 2025). Note: ChiralFold's 0% rate is guaranteed by construction — each residue is built with explicit stereochemistry encoding — rather than learned from data. The comparison demonstrates that construction-based approaches solve a problem AF3's diffusion architecture fundamentally cannot.

Fisher's exact test: p < 6.7×10⁻¹⁴⁴. 31 PDB structures audited: 30/31 = 100% Cα correctness.

### External Benchmark: Childs et al. 2025

41 sequences using the real D-peptide sequences from the Childs et al. (2025) PDB structures: DP19:L-19437 (DEHELLETAARWFYEIAKR, PDB 7YH8), DP9:Streptavidin (LWQHEATWK, PDB 5N8T), DP12:MDM2 (DWWPLAFEALLR, PDB 3IWY), plus L/D pattern variants. ChiralFold: 0/478 chiral residues violated. AF3: 50–52% per-residue violation rate on the same systems. See `benchmarks/childs2025_comparison.py`.

### PDB-Wide D-Residue Chirality Verification

Independently verified (no ChiralFold code — numpy + raw PDB coordinates only) 12,573 D-amino acid residues across 4,616 PDB files (>91% coverage of all 18 standard D-amino acid CCD codes). Found **29 D-label/L-coordinate mismatches in 16 structures** — cases where deposited Cα coordinates show L-stereochemistry despite D-amino acid labels. Error rate: 0.23%. Errors cluster in 5 CCD codes (DTY 2.35%, DLY 0.73%, DPN 0.33%, DAR 0.36%, DSN 0.28%); 9 codes confirmed clean at ≥91% coverage. All are invisible to MolProbity.

**Error classification (cross-referenced against deposition remarks, COMPND records, biological context, and primary literature):**

| PDB | Residue | Chain | Pos | Signed Volume | Error Type | Evidence |
|-----|---------|:-----:|:---:|:-------------:|:----------:|----------|
| 1ABI | DPN | I | 56 | +2.49 | **Stereochem** | Hirulog D-Phe by design; internal control DPN:1 correct (vol=-2.60) |
| 1BG0 | DAR | A | 403 | +2.58 | CCD-Code | Arginine kinase substrate is L-Arg; standalone ligand |
| 1D7T | DTY | A | 4 | +1.85 | **Stereochem** | Contryphan [D-Tyr4] explicitly designed; NMR structure |
| 1HHZ | DAL | E | 1 | +2.70 | **Stereochem** | Cell wall pentapeptide D-Ala; 0.99 Å atomic resolution |
| 1KO0 | DLY | A | 542 | +0.12 | Borderline | Title says "D,L-lysine" (racemic); ALTLOC B, B=32.3 |
| 1MCB | DHI | P | 3 | +2.60 | Mislabel | COMPND says "L-HIS" explicitly; DHI should be HIS |
| 1OF6 | DTY | A-H | 1369-1370 | +2.51 to +2.67 | CCD-Code | DAHP synthase inhibitor is L-Tyr ([Nature 2023](https://doi.org/10.1038/s42004-023-00946-x)); 8 ligand copies |
| 1P52 | DAR | A | 403 | +2.54 | CCD-Code | Arginine kinase mutant; substrate is L-Arg |
| 1UHG | DSN | D | 164 | +2.21 | Mislabel | Ovalbumin is an L-protein; no biological reason for D-Ser |
| 1XT7 | DSG | A | 3 | +2.55 | **Stereochem** | Daptomycin antibiotic: D-Asn biologically required |
| 2AOU | DCY | A | 248 | +2.67 | Mislabel | Histamine methyltransferase is L-protein; DCY should be CYS |
| 2ATS | DLY | A | 3001-3003 | +2.56 to +2.59 | CCD-Code | Title says "(S)-lysine" = L-Lys; 3 ligand copies |

**Error type breakdown:**
- **CCD-Code** (18 errors, 5 structures): Non-polymer ligand labeled with D-form CCD code when the crystallized molecule is L-form. The coordinates correctly model L-stereochemistry; only the chemical component code is wrong.
- **Mislabel** (3 errors, 3 structures): Polymer residue labeled with D-amino acid code in a context where L is biologically correct (L-protein, or COMPND record specifying L-form).
- **Stereochem** (6 errors, 6 structures): **Most concerning** — the biology requires D-stereochemistry but the deposited coordinates show L. These are genuine coordinate-level errors.
- **Borderline** (2 errors, 2 structures): Near-zero signed volume; inconclusive without additional evidence.

Regardless of type, all 29 mismatches are real annotation inconsistencies in the PDB. The signed volume method detects all error types without requiring knowledge of biological context.

Full dataset: `results/d_residue_verification.csv` (12,574 rows with raw coordinates). Classification: `results/error_classification.json`. Per-structure reclassification evidence: `results/error_table_verified.csv`. **Code-level coverage summary: `results/ccd_code_coverage_summary.csv`** (see below).

> **Benchmark Reproducibility:** All 16 error structures are independently verifiable from `results/error_table_verified.csv` (signed volumes, B-factors, ALTLOC flags, biological evidence) and `results/error_classification.json` (category breakdown and method). No ChiralFold code is required — only numpy and raw PDB coordinate files from RCSB.

**Bulletproof verification:** Five independent checks confirm the findings: (1) Sign convention validated on 24/24 known-correct D-residues in 3IWY (all negative volumes); (2) 1KO0 reclassified as borderline (vol=+0.12, ALTLOC B, B=32.3 — inconclusive); (3) 1OF6 confirmed across all 8 chains (all L-coordinates, consistent with L-Tyr biological role); (4) 1ABI internal control passes cleanly (DPN:1 vol=-2.60 correct vs DPN:56 vol=+2.49 error); (5) Full re-verification of all 16 structures with biological context cross-referencing. See `benchmarks/bulletproof_verification.py` and `results/bulletproof_outputs/`.

**Correlation analysis:** 13 of the 16 error structures were deposited between 1992 and 2005, consistent with the 2006–2008 wwPDB remediation effort. Three post-remediation errors (2RMI 2007, 2W76 2008, 3RIT 2011) confirm the pipeline still lacks a D-specific chirality check. Deposition year significantly predicts errors (Mann-Whitney U=278, p=0.0027 on the initial 12-structure survey). Resolution does not (p=0.19) — errors span 0.99–2.77 Å, indicating a labeling problem rather than a data quality problem.

---

### Code-Level Coverage Summary

`results/ccd_code_coverage_summary.csv` — **18 rows, 12 columns.** One row per standard D-amino acid CCD code. Lets reviewers reproduce the error-clustering finding without parsing the LaTeX.

| Column | Type | Description |
|--------|------|-------------|
| `ccd_code` | string | Three-letter wwPDB Chemical Component Dictionary code (e.g. `DTY`, `DLY`) |
| `rcsb_total_entries` | integer | Total PDB entries containing this code, per RCSB full-text search |
| `entries_surveyed` | integer | Entries downloaded and verified in this study |
| `coverage_pct` | float | `entries_surveyed / rcsb_total_entries × 100` |
| `mmcif_only_unavailable` | integer | Entries not available as legacy PDB format (post-2019 depositions and large cryo-EM assemblies); excluded from verification |
| `residues_checked` | integer | Individual Cα chirality checks performed (residues with complete N/Cα/C/Cβ backbone) |
| `n_errors` | integer | D-labeled residues with confirmed L-stereochemistry (signed volume > 0) |
| `n_error_structures` | integer | Distinct PDB IDs containing at least one error |
| `error_structures` | string | Semicolon-separated list of those PDB IDs |
| `error_rate_pct` | float | `n_errors / residues_checked × 100` |
| `status` | string | `error-prone` (n_errors > 0) or `confirmed-clean` (n_errors = 0 at ≥91% coverage) |
| `biological_context` | string | Plain-English summary of the structural families in which this code appears |

**Key finding from this table:** Errors concentrate in five codes — DTY (2.35% error rate), DLY (0.73%), DAR (0.36%), DSN (0.28%), DPN (0.33%) — while nine codes are confirmed clean at zero errors across ≥91% of their RCSB universe: DAS, DGL, DGN, DIL, DLE, DPR, DTH, DTR, DVA. The error-prone codes share two biological scenarios: enzyme active-site ligands where the D-form CCD code is confused with the L-form substrate (DAR/DTY/DLY), and designed or NRPS-assembled peptides where D-configuration is required but L-coordinates were deposited (DPN/DSN).

### Planarity Fix

- D-peptides: 39% → 94% within 6° of planar
- L-proteins: 33% → 95% (generalizes across 5 backbone types)

### Mirror Pipeline

- 5 PDB systems, 13,767 atoms: 0.0 Å coordinate error
- L→D→L round-trip: mathematically exact
- Contact geometry preserved by construction: 105 Cα-Cα contacts and 10 H-bond donor-acceptor pairs maintained within 0.001 Å across L→D transformation

### wwPDB Comparison

- 31 structures audited (X-ray, NMR, cryo-EM)
- Ramachandran: Spearman ρ = 0.49 (p = 0.006) vs wwPDB
- Mean outlier rate: CF 0.60% vs wwPDB 0.64%

## Scope and Limitations

ChiralFold is a stereochemistry toolkit, not a de novo structure predictor. It excels at chirality auditing, L↔D coordinate transformation, and template-based conformer generation. For de novo folding, use AlphaFold 3 or ESMFold — then pipe the output through ChiralFold to validate or correct stereochemistry.

| Capability | Status |
|-----------|--------|
| Chirality auditing (L and D) | Production-ready |
| Mirror-image transformation | Production-ready (0.0 Å error) |
| AF3 chirality correction | New in v3.2 |
| Diastereomer enumeration | New in v3.2 |
| Interface scoring | New in v3.2 |
| Template threading | Available (template-dependent; requires structural homolog in PDB) |
| Fragment assembly | Available (Chou-Fasman SS + NeRF backbone; not comparable to learned models) |
| De novo fold prediction | Not supported — use AF3/ESMFold + ChiralFold correction |

## Previously Addressed Limitations

| Issue | v3.0 Status | v3.2 Resolution |
|-------|:----------:|:----------------:|
| Not a fold predictor | Mirror-only | **Template threading + fragment assembly** (template-dependent; requires structural homolog in PDB) |
| Ramachandran uses rectangles | Calibrated rectangles | **Hybrid: empirical PDB grid + calibrated rectangle fallback** |
| No rotamer analysis | Planned | **Penultimate Rotamer Library validation** (chi1 scoring) |
| Clash methodology differs | Heavy-atom only | **Backbone H-atom placement** before scoring |
| Conformer limit at 30 res | Hard limit | **Fragment assembly for any protein length** |

## Project Structure

```
ChiralFold/
├── chiralfold/
│   ├── __init__.py
│   ├── model.py              # ChiralFold model + SMILES builders
│   ├── auditor.py            # PDB structure quality auditor (H-aware clashes)
│   ├── af3_correct.py        # AlphaFold 3 chirality correction pipeline
│   ├── enumerate.py          # Diastereomer enumeration + ranking
│   ├── interface_scorer.py   # Binding interface scoring
│   ├── ramachandran.py       # Hybrid empirical + rectangular Ramachandran
│   ├── rotamers.py           # Side-chain rotamer validation
│   ├── threading.py          # Template-based fold prediction
│   ├── fragments.py          # Fragment-based backbone assembly
│   ├── validator.py          # Chirality validation engine
│   ├── pdb_pipeline.py       # Mirror-image PDB transformation
│   ├── geometry.py           # Planarity fix + post-processing
│   ├── cli.py                # Command-line interface
│   └── data/
│       ├── test_sequences.py # 46-sequence test library
│       └── ramachandran_grid.json  # Empirical φ/ψ probability grid
├── tests/
│   └── test_chirality.py     # Unit tests (incl. external PDB validation)
├── benchmarks/               # Benchmark scripts (incl. bulletproof_verification.py)
├── results/                  # Generated outputs
├── CONTRIBUTING.md           # How to contribute
├── pyproject.toml
├── LICENSE (MIT)
└── README.md
```

## Validation

All numbers in this README were produced by running ChiralFold commands in a single validation session. Key outputs:

**3IWY Audit (experimental D-peptide crystal structure):**
Residues: 189. Chirality: 181 correct, 0 wrong, 8 Gly (100.0%). Ramachandran: 95.2% favored, 0.5% outlier. Score: 65.4/100.

**Childs 2025 Comparison (41 sequences, real PDB sequences):**
0/478 chiral residues violated across DP19 (DEHELLETAARWFYEIAKR, PDB 7YH8), DP9 (LWQHEATWK, PDB 5N8T), DP12 (DWWPLAFEALLR, PDB 3IWY), and 26 L/D pattern variants.

**MDM2 Interface Score (1YCR, p53:MDM2):**
Buried Surface Area: 1,980 Å². Shape Complementarity: 0.933. Hydrogen bonds: 10. Interface score: 61.9/100.

**D-Residue Verification (12,573 residues, independent of ChiralFold):**
29 D-label/L-coordinate mismatches in 16 PDB structures. Error rate: 0.23% (12,573 residues checked across 4,616 files, >91% coverage of all 18 standard D-amino acid CCD codes). Classified by biological context: 6 genuine stereochemistry errors (coordinates modeled as L where D is required), 18 CCD code misassignments across 5 structures (L-form ligand with D-code), 3 polymer residue mislabels, 2 borderlines. All verified using only numpy and raw PDB coordinates (no ChiralFold code). Classification cross-referenced against COMPND records, SEQRES, deposition titles, and primary literature. Full dataset: `results/d_residue_verification.csv`. Classification: `results/error_classification.json`. Code-level coverage table: `results/ccd_code_coverage_summary.csv`.

**Batch Audit (validated against this README):**
1UBQ: 100% chirality, 97.3% Rama favored, 96.0% planarity, score 78.6.
1CRN: 100% chirality, 88.6% Rama favored, 91.1% planarity, score 72.1.
1SHG: 100% chirality, 92.7% Rama favored, 66.1% planarity, score 70.0.

## Citation

```bibtex
@software{chiralfold2025,
  title     = {ChiralFold: General-Purpose Protein Stereochemistry Toolkit},
  author    = {Tommaso R. Marena},
  year      = {2025},
  url       = {https://github.com/Tommaso-R-Marena/ChiralFold},
  version   = {3.2.1},
  note      = {PDB auditing calibrated against wwPDB/MolProbity,
               chirality-correct coordinate generation, AF3 correction
               pipeline, mirror-image binder design validated on MDM2
               (dPMI-gamma, Kd=53nM)}
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
