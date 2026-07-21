# Bioinformatics submission highlights (for editorial form / cover letter)

1. First systematic audit of D-amino acid stereochemistry across 12,573 PDB residues (>91% RCSB coverage per CCD code).

2. 29 D-label/L-coordinate mismatches in 16 structures escape MolProbity's L-only Cα chirality check; six are genuine stereochemistry errors.

3. ChiralFold auditor agrees with wwPDB Ramachandran outliers on 362 era-representative structures (Spearman rho = 0.52).

4. Machine-checked Lean 4 proofs: chirality sign cannot factor through pairwise distances alone.

5. Open software: pip install chiralfold, Docker, GitHub v3.5.1, HF Space (Light/Dark; RCSB + AlphaFold DB fetch; PDB/mmCIF/FASTA), full reproduction package.

- Full live mmCIF-only D-residue universe survey (79 entries): recovers 29/29 known errors; finds 9BC4 DLE mismatch.

- Machine-checked Lean 4 chirality no-go proofs with full derivations in Supplementary Notes S1.

- Clashscore uses topology-aware 1–2/1–3/1–4 exclusions (secondary UI metric; survey/Rama claims unaffected; panel mean cf_clash 23.8 vs wwPDB 18.0).
