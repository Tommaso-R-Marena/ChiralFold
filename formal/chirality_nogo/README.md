# Formal Chirality No-Go Theorems (Lean 4)

Machine-checked proofs that **pairwise distances alone cannot determine chirality sign**, and that any model whose output factors through the pairwise distance matrix is necessarily reflection-invariant. These results provide a formal foundation for ChiralFold's design choice to use signed tetrahedron volumes at Cα rather than distance-only geometry.

Formalized with [Aristotle](https://aristotle.harmonic.fun) (Harmonic) on Lean 4 v4.28.0 + Mathlib.

## Main results

| File | Key theorem | Statement |
|------|-------------|-----------|
| `ChiralityNoGo.lean` | `same_distances_opposite_orientations_standardTet` | A tetrahedron and its mirror share identical pairwise squared distances but have opposite oriented volumes |
| `ChiralityNoGo.lean` | `no_distance_only_classifier_separates_standardTet` | No Boolean classifier on distance matrices can label a configuration and its mirror differently |
| `OrthogonalChirality.lean` | `improper_orthogonal_distance_chirality_obstruction` | Improper orthogonal maps (det = −1) preserve all distances but negate oriented volume |
| `ModelAbstraction.lean` | `orientTarget_not_factorsThroughDistMatrix` | The oriented volume target cannot be represented by any distance-factoring model |
| `ModelAbstraction.lean` | `reflection_odd_target_not_distance_factorable` | Any target that flips sign under reflection cannot factor through distances |
| `DistancesOrientVolume.lean` | `orient_sq_determined_by_pairwise_distances` | Squared oriented volume is determined by pairwise distances (via Gram matrix determinant) |
| `DistancesOrientVolume.lean` | `distance_data_determines_magnitude_not_chirality_sign` | Same distances imply equal \|orient\|² but not equal orient sign |

All proofs compile without `sorry` or `admit`, using only standard Lean axioms (`propext`, `Classical.choice`, `Quot.sound`).

## Build

Requires [Lean 4 v4.28.0](https://leanprover.github.io/lean4/doc/setup.html) and [Lake](https://github.com/leanprover/lean4/tree/master/src/lake).

```bash
cd formal/chirality_nogo
lake exe cache get   # optional: fetch prebuilt Mathlib cache
lake build
```

## Relation to ChiralFold

ChiralFold's Cα chirality check (Eq. 1 in the paper) computes a **signed oriented volume** from four backbone atoms (N, Cα, C, Cβ) — equivalent to the `orient4` construction in `ChiralityNoGo.lean`. This is strictly stronger than pairwise-distance information: the formal proofs show that distance-only models (including EGNN-style architectures and diffusion denoisers that operate on interatomic distances) cannot, in principle, distinguish enantiomers without additional signed-geometry or stereochemical inputs. ChiralFold's construction-based pipeline (explicit SMILES stereochemistry, signed-volume auditing, AF3 reflection correction) is designed around this obstruction.

## Citation

When referencing these proofs, cite the ChiralFold paper and acknowledge Aristotle:

```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

See `ARISTOTLE_SUMMARY.md` for per-run theorem inventories.
