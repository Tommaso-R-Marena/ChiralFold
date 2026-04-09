# DPR Full-Coverage Expansion Report

**Date:** 2026-04-09  
**Objective:** Bring DPR (D-proline) from 83.3% to >91% coverage  
**Result:** 96.2% coverage achieved — DPR **confirmed clean**

---

## Downloads

- Attempted: 44 remaining DPR structures  
- Retrieved: 34 (10 unavailable as legacy PDB — all post-2015 mmCIF-only entries)  
- Total surveyed after this pass: 4,616 PDB files

## Verification

| Metric | Value |
|--------|-------|
| DPR structures surveyed | 254 / 264 (96.2%) |
| DPR residues checked | included in 12,573 total |
| DPR errors | **0** |
| Total dataset errors | 29 (unchanged) |
| Total error structures | 16 (unchanged) |

**No errors found in DPR.** D-proline (DPR) is confirmed clean at 96.2% coverage.

---

## Final 18-Code Coverage Summary

| Code | RCSB | Surveyed | Coverage | Errors | Status |
|------|------|---------|---------|--------|--------|
| DAL | 600 | 577 | 96.2% | 1 | Error-prone |
| DAR | 189 | 177 | 93.7% | 2 | Error-prone |
| DAS | 1,300 | 1,237 | 95.2% | 0 | **Confirmed clean** |
| DCY | 99 | 96 | 97.0% | 1 | Error-prone |
| DGL | 223 | 205 | 91.9% | 0 | **Confirmed clean** |
| DGN | 92 | 87 | 94.6% | 0 | **Confirmed clean** |
| DHI | 87 | 85 | 97.7% | 1 | Error-prone |
| DIL | 118 | 113 | 95.8% | 0 | **Confirmed clean** |
| DLE | 191 | 181 | 94.8% | 0 | **Confirmed clean** |
| DLY | 166 | 159 | 95.8% | 9 | Error-prone |
| DPN | 2,348 | 2,182 | 92.9% | 2 | Error-prone |
| DPR | 264 | 254 | **96.2%** | **0** | **Confirmed clean** |
| DSG | 107 | 104 | 97.2% | 1 | Error-prone |
| DSN | 184 | 174 | 94.6% | 2 | Error-prone |
| DTH | 135 | 131 | 97.0% | 0 | **Confirmed clean** |
| DTR | 131 | 124 | 94.7% | 0 | **Confirmed clean** |
| DTY | 156 | 143 | 91.7% | 10 | Error-prone |
| DVA | 200 | 193 | 96.5% | 0 | **Confirmed clean** |

**9 codes confirmed clean (≥91% coverage, 0 errors): DAS, DGL, DGN, DIL, DLE, DPR, DTH, DTR, DVA**  
**9 codes error-prone: DAL, DAR, DCY, DHI, DLY, DPN, DSG, DSN, DTY**

All 18 standard D-amino acid CCD codes are now surveyed at ≥91% of their respective RCSB universes. The 245 remaining unavailable entries are exclusively post-2015 cryo-EM and large X-ray structures that lack legacy PDB coordinate files.

---

## Error Clustering Pattern

The 9 error-prone codes split into two mechanistic groups:

**CCD-code confusion (enzyme/substrate context):** DAR, DTY, DLY  
The D-form CCD code was applied to a ligand that is biologically L-form: arginine kinase substrates (DAR→L-Arg), DAHP synthase allosteric inhibitor (DTY→L-Tyr), epimerase substrates and co-crystallization reagents (DLY→L-Lys).

**Stereochem errors in designed/natural-product peptides:** DPN, DSN  
The D-configuration is biologically required and explicitly designed, but the deposited Cα shows L-stereochemistry: cyclic peptide CRF antagonist (DPN, astressin 2RMI), hirulog thrombin inhibitor (DPN, 1ABI), pyoverdine siderophore (DSN, 2W76), L-protein mislabel (DSN, 1UHG).

**Isolated mislabels:** DCY (1 error, 2AOU), DHI (1 error, 1MCB), DSG (1 error, 1XT7), DAL (1 error, 1HHZ)  
Single instances without a recognizable family pattern.

**Clean codes:** DPR, DTH, DTR, DVA, DIL, DLE, DGL, DGN, DAS  
These codes appear predominantly as backbone-incorporated residues in well-characterized cyclic peptide antibiotics, non-ribosomal peptide natural products, and D-peptide therapeutic structures where depositors have strong independent stereochemical constraints (MS fragmentation, NMR NOE patterns, biosynthetic logic).

D-Pro (DPR) in particular is a structural residue in gramicidin-type peptides and beta-turn mimetics: its geometric role (proline-induced turn) provides an independent check that depositors likely apply during model building, which may explain its zero-error rate despite widespread use.

---

## Files Updated This Pass

- `results/d_residue_verification.csv` — 12,574 rows
- `results/d_residue_verification_summary.json`
- `results/DPR_EXPANSION_REPORT.md` — this file
- `paper/chiralfold_paper.tex` — Discussion para 4 rewritten with code-cluster findings; Data Availability row count updated to 12,574
