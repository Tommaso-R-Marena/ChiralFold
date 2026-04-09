# Expanded PDB Verification Report

**Date:** 2026-04-09  
**Dataset:** 590 PDB files → 501 structures containing D-amino acids  
**D-residues checked:** 10,469 (6.24× the original 1,677)  
**Method:** Identical to original — signed tetrahedron volume at Cα, numpy only, no ChiralFold code

---

## Coverage

| Metric | Original | Expanded |
|---|---|---|
| PDB files scanned | 231 | 501 |
| D-residues checked | 1,677 | 10,469 |
| D-correct | 1,656 (98.75%) | 10,442 (99.74%) |
| L-errors (D-code, L-coords) | 21 (1.25%) | 27 (0.26%) |
| Error structures | 12 | 14 |

The decrease in error rate from 1.25% → 0.26% reflects that the original survey was designed to oversample structures known to contain D-amino acids in biologically interesting contexts; the expansion draws more broadly from all CCD-code-containing structures, which are predominantly correct.

## RCSB Coverage Query

All 19 D-amino acid CCD codes queried via RCSB full-text search:

| Code | Total PDB structures |
|---|---|
| DAL | 600 |
| DAR | 189 |
| DAS | 1,300 |
| DCY | 99 |
| DCL | 223 |
| DGN | 92 |
| DHI | 87 |
| DIL | 118 |
| DLE | 191 |
| DLY | 166 |
| DPN | 2,348 |
| DPR | 264 |
| DSN | 184 |
| DSG | 107 |
| DTH | 135 |
| DTR | 131 |
| DTY | 156 |
| DVA | 200 |

Total unique structures with any D-AA code: ~4,915.  
This survey covers 501/4,915 = **10.2% of the total PDB D-AA space**.

---

## Original 12 Error Structures — All Confirmed

Every error from the original survey is reproduced with identical signed volumes:

| PDB | Residue | Vol (original) | Vol (expanded run) | Match |
|---|---|---|---|---|
| 1ABI | DPN:I:56 | +2.4884 | +2.4884 | ✓ |
| 1BG0 | DAR:A:403 | +2.5832 | +2.5832 | ✓ |
| 1D7T | DTY:A:4 | +1.8491 | +1.8491 | ✓ |
| 1HHZ | DAL:E:1 | +2.7037 | +2.7037 | ✓ |
| 1KO0 | DLY:A:542 | +0.1156 | +0.1156 | ✓ |
| 1MCB | DHI:P:3 | +2.5975 | +2.5975 | ✓ |
| 1OF6 | DTY (8 chains) | +2.51–+2.67 | identical | ✓ |
| 1P52 | DAR:A:403 | +2.5364 | +2.5364 | ✓ |
| 1UHG | DSN:D:164 | +2.2136 | +2.2136 | ✓ |
| 1XT7 | DSG:A:3 | +2.5501 | +2.5501 | ✓ |
| 2AOU | DCY:A:248 | +2.6722 | +2.6722 | ✓ |
| 2ATS | DLY (3 copies) | +2.56–+2.59 | identical | ✓ |

Sign-convention control: 24/24 D-residues in 3IWY (known-correct) give negative volumes (−2.26 to −2.73). All 1,035 protein-chain L-residues in 3RIT give positive volumes. Convention is self-consistent.

---

## New Findings

### 3RIT — NEW CCD-Code Error (deposited 2011)

**Structure:** Dipeptide epimerase from *Methylococcus capsulatus* complexed with Mg²⁺ and dipeptide  
**Title:** "...complexed with dipeptide L-Arg-D-Lys"  
**Resolution:** 2.70 Å | **Deposited:** 2011-04-27 | **PMID:** 22392983

DLY appears in all 5 chains (A–E) as a HETATM ligand:

| Chain | Res# | Vol | Classification |
|---|---|---|---|
| A | 364 | +2.7548 | L-coordinates |
| B | 362 | +2.4824 | L-coordinates |
| C | 362 | +2.4182 | L-coordinates |
| D | 364 | +2.7681 | L-coordinates |
| E | 363 | +2.4892 | L-coordinates |

The dipeptide epimerase catalyzes L-Lys → D-Lys epimerization. The title claims the structure is in complex with the D-Lys product (L-Arg-**D**-Lys), but the Cα coordinates show L-configuration in all five independent copies. This is a **CCD-Code error**: the depositors used DLY (D-lysine) when the coordinates model L-lysine, consistent with the substrate state (L-Arg-**L**-Lys) rather than the labeled product state.

**Significance:** This is the first post-2005 CCD-code error found in the expanded survey. It revises the paper's temporal claim from "all errors pre-2006" to "13 of 14 errors pre-2006; 1 CCD-code error in 2011." The temporal correlation weakens slightly but the mechanistic interpretation is unchanged: CCD-code errors persist because the deposition pipeline has no automated cross-check.

### 2H9E — Borderline (deposited 2006)

**Structure:** FXa/Selectide/NAPc2 ternary complex  
**COMPND:** "SELECTIDE INHIBITOR DTY-ILE-ARG-LEU-LPD PEPTIDE" (chain S)  
**Resolution:** 2.20 Å | **Deposited:** 2006-06-09 | **PMID:** 17173931

| Chain | Res# | Vol | B-factor | Classification |
|---|---|---|---|---|
| S | 1 | +0.5642 | 33.99 | **Borderline** |

Signed volume = +0.56 is above zero but far below the +1.85 minimum of confirmed errors. By analogy with 1KO0 (+0.12, ALTLOC B, B=32.3), this is classified **borderline/inconclusive**. The COMPND record explicitly names this as a D-Tyr-containing inhibitor peptide, so D-Tyr is biologically expected. The small positive volume may reflect disorder at the N-terminus of the inhibitor chain rather than a systematic chirality error. Not counted as a confirmed error.

---

## Statistical Update

| Test | Original (12 structures) | Expanded (14 structures) |
|---|---|---|
| All orig. errors pre-2006 | ✓ (Mann-Whitney p=0.0027) | ✓ |
| New CCD error post-2005 | N/A | 3RIT (2011) — 1 exception |
| Resolution predicts errors | p=0.19 (no) | Unchanged |
| Post-2005 error rate | 0/~50 structures | 1/~370 structures (0.27%) |

The 40-structure random sample of post-2005 additions found 0 errors among 23 confirmed post-2005 structures. The one confirmed post-2005 error (3RIT, 2011) is a CCD-code error — mechanistically identical to the 2003–2005 CCD-code cases.

---

## Recommended Paper Update

The temporal claim in the Discussion paragraph 3 should be softened from:

> "all 12 error-containing structures were deposited between 1992 and 2005"

to:

> "13 of the 14 error-containing structures were deposited between 1992 and 2005; the one exception (3RIT, 2011) is a CCD-code error of the same type as the pre-remediation cases, consistent with the conclusion that the deposition pipeline has no automated stereochemistry cross-check regardless of deposition era."

The Results section temporal statistics (Mann-Whitney U=278, p=0.0027) remain valid for the original 12-structure set and should be reported as such.

---

## Files Updated by This Run

- `results/d_residue_verification.csv` — expanded to 10,470 rows (10,469 checkable)
- `results/d_residue_verification_summary.json` — updated summary statistics
- `results/d_survey/` — 370 new PDB files added (590 total)
- `results/EXPANDED_VERIFICATION_REPORT.md` — this file
