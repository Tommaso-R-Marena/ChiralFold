# Bioinformatics submission highlights (for editorial form / cover letter)

1. First systematic audit of D-amino acid stereochemistry across 12,573 PDB residues (>91% RCSB coverage per CCD code).

2. 29 D-label/L-coordinate mismatches in 16 structures escape MolProbity's L-only Cα chirality check; six are genuine stereochemistry errors.

3. ChiralFold auditor agrees with wwPDB Ramachandran outliers on 362 era-representative structures (Spearman rho = 0.52).

4. Machine-checked Lean 4 proofs: chirality sign cannot factor through pairwise distances alone.

5. Deployable at archive scale: the auditor runs at 15-18 us per atom with wall time linear in structure size (1,589-25,424 atoms), 3.5-3.9x faster than the previous release; interface scoring 3.1-6.1x faster with peak memory cut from 1,260 MB to 183 MB. Offline harnesses reproduce both the absolute and before/after numbers.

6. Open software: pip install chiralfold, Docker, GitHub v3.6.0, HF Space (Light/Dark; RCSB + AlphaFold DB fetch; PDB/mmCIF/FASTA), full reproduction package.

- Full live mmCIF-only D-residue universe survey (79 entries): recovers 29/29 known errors; finds 9BC4 DLE mismatch.

- Machine-checked Lean 4 chirality no-go proofs with full derivations in Supplementary Notes S1.

- Clashscore: topology-aware exclusions + Probe-style HA/HA2/HA3 + geometry-gated H-bonds (secondary UI metric; survey/Rama unaffected; panel mean cf_clash 27.4 vs wwPDB 18.0).
