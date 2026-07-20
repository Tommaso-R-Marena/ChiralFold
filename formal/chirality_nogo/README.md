This project was edited by [Aristotle](https://aristotle.harmonic.fun).

To cite Aristotle:
- Tag @Aristotle-Harmonic on GitHub PRs/issues
- Add as co-author to commits:
```
Co-authored-by: Aristotle (Harmonic) <aristotle-harmonic@harmonic.fun>
```

# Formal Chirality No-Go Theorems (Lean 4)

This project machine-checks a precisely scoped obstruction for **distance-only**
geometric representations. For four points in `ℝ³`, the six pairwise squared
distances determine the square (and hence magnitude) of the oriented tetrahedral
volume, but not its sign. Consequently, a model or classifier whose entire input
factors through that distance matrix cannot distinguish a configuration from its
mirror image.

## What is proved

- Every orthogonal transformation preserves the pairwise distance matrix.
- Oriented volume is multiplied by the determinant of the orthogonal map.
  Thus every improper orthogonal map (`det = -1`) reverses signed orientation.
- Every output that factors through `distMatrix4` is invariant under orthogonal
  transformations, including improper ones.
- `orient4` and `chiralitySignBool` do not factor through `distMatrix4`.
- A continuous real-vector representation that factors through distances cannot,
  after any downstream classifier, agree with chirality sign on all configurations.
- For atoms ordered `(N, Cα, C, Cβ)`, their signed tetrahedral volume cannot be a
  function of the six pairwise distances alone.
- The Gram determinant equals squared oriented volume, so distances recover
  magnitude but not sign.

## What is **not** proved

This is not an impossibility theorem for AlphaFold 3 or for all SE(3)-equivariant
neural networks. An SE(3)-equivariant architecture may use angular, tensorial,
higher-order, sequence-derived, or other orientation-sensitive information and
need not factor through a pairwise distance matrix. The formal bridge theorem
therefore states its full-orthogonal-invariance hypothesis explicitly. No claim is
made that arbitrary GNNs or diffusion predictors satisfy that hypothesis.

The result supports signed-volume auditing as a minimal way to expose information
that pairwise distances omit; it does not explain or prove a failure of a complete
predictor architecture.

## Suggested manuscript sentence

> Machine-checked Lean 4 proofs show that the six pairwise distances among four
> atoms determine the magnitude but not the sign of their signed tetrahedral
> volume, so representations factoring solely through those distances are blind
> to enantiomeric configuration; this distance-only obstruction is not an
> impossibility result for SE(3)-equivariant predictors with angular or
> higher-order features, but motivates explicit signed-volume auditing.

## Main declarations

| File | Declaration | Formal consequence |
|---|---|---|
| `ChiralityNoGo.lean` | `same_distances_opposite_orientations_standardTet` | A tetrahedron and its mirror have equal distances and opposite orientation. |
| `OrthogonalChirality.lean` | `improper_orthogonal_distance_chirality_obstruction` | Any orthogonal map with determinant `-1` preserves distances and flips orientation. |
| `ModelAbstraction.lean` | `factorsThroughDistMatrix_orthogonal_invariant` | Every distance-factoring output is invariant under every orthogonal map. |
| `ModelAbstraction.lean` | `orient4_not_distance_factorable_via_improper_orthogonal` | An arbitrary improper orthogonal map supplies the signed-orientation obstruction at every nondegenerate witness. |
| `ModelAbstraction.lean` | `chiralitySignBool_not_distance_factorable` | Boolean chirality sign cannot be recovered from the distance matrix. |
| `ModelAbstraction.lean` | `no_continuous_distance_representation_classifier` | Continuity does not rescue a representation that still factors through distances. |
| `ModelAbstraction.lean` | `proteinCaOrient_not_function_of_six_pairwise_distances` | `(N, Cα, C, Cβ)` signed volume is not a function of the six distances alone. |
| `DistancesOrientVolume.lean` | `det_gram3_eq_orient_sq` | The Gram determinant equals squared oriented volume. |
| `Examples.lean` | `standardTetExample_same_distances_different_orientations` | A concrete normalized example exhibits the ambiguity. |

Existing theorem names used by earlier supplementary material are preserved.

## Build

```bash
lake build
```

The production Lean files contain no `sorry` or `admit`.

## Paper-facing text

A ready-to-paste paragraph and theorem-to-biological-consequence table are in
[`paper_methods_snippet.tex`](paper_methods_snippet.tex).
