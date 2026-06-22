# Summary of changes for run 79369806-0d7e-4b09-9aaa-ef6fc575bb60
Created `RequestProject/ModelAbstraction.lean` with all requested theorems, fully proven (no `sorry`, no `admit`, no non-standard axioms).

## Theorem names

1. **`factorsThroughDistMatrix_reflection_invariant`** — Any model `M : Config4 → α` that factors through the pairwise distance matrix satisfies `M (reflectZConfig P) = M P` for all configurations `P`. Proof: unpack the factoring witness `F`, rewrite both sides using `distMatrix4_reflectZConfig`.

2. **`orientTarget_not_factorsThroughDistMatrix`** — The signed orientation target `orientTarget P = orient4 P` cannot factor through the distance matrix. Proof: if it did, reflection invariance would give `orient4(reflect standardTet) = orient4(standardTet)`, but one is `1` and the other is `-1`.

3. **`reflection_odd_target_not_distance_factorable`** — General theorem: any target `T : Config4 → ℝ` satisfying `T(reflect P) = -T(P)` and `T(P) ≠ 0` for some `P` cannot factor through distances. Proof: factoring implies `T(P) = T(reflect P) = -T(P)`, so `2·T(P) = 0`, hence `T(P) = 0`, contradicting `T(P) ≠ 0`.

4. **`orient4_not_distance_factorable_general`** — Instantiation of the general theorem to `orient4` using `standardTet` as witness.

5. **`chiralitySignBool_not_distance_factorable`** — The Boolean chirality classifier `chiralitySignBool` (true iff `orient4 P > 0`) cannot factor through distances. Proof: it evaluates to `true` on `standardTet` and `false` on its reflection, but distances are identical.

## Key definition

- **`FactorsThroughDistMatrix M`**: `∃ F, ∀ P, M P = F (distMatrix4 P)` — a model's output depends only on pairwise distances.

## Summary (under 250 words)

The file formalizes a no-go theorem for distance-based ML models applied to chirality detection. The central definition, `FactorsThroughDistMatrix`, captures any model whose output is a function of the pairwise distance matrix. The main result shows that such models are necessarily reflection-invariant (`factorsThroughDistMatrix_reflection_invariant`), because reflecting a configuration preserves all pairwise distances. Consequently, any target function that changes sign under reflection — such as the oriented volume `orient4` — cannot be represented by any distance-factoring model. This is proven both for the specific `orientTarget` function and in full generality via `reflection_odd_target_not_distance_factorable`. A Boolean version (`chiralitySignBool_not_distance_factorable`) demonstrates that even a binary chirality classifier cannot factor through distances. All proofs compile without sorry and use only standard axioms (propext, Classical.choice, Quot.sound).

# Summary of changes for run c4ae68de-fa88-4ee9-95fc-41267543c062
Created `RequestProject/DistancesOrientVolume.lean` — a complete, sorry-free Lean 4 file proving that pairwise squared distances determine squared oriented volume but not its sign. All theorems compile and depend only on standard axioms.

**Theorems proved:**

1. **`distSq_eq_normSq_vsub`** — distSq equals the squared norm of the difference vector.
2. **`gram3_eq_transpose_mul_edgeMatrix`** — the Gram matrix G = Mᵀ·M where M is the edge matrix.
3. **`orient_eq_det_edgeMatrixFromPoints`** — orient(a,b,c,d) = det(edge matrix).
4. **`det_gram3_eq_orient_sq`** — det(G) = orient(a,b,c,d)², via det(MᵀM) = det(M)².
5. **`dot_vsub_vsub_eq_distSq_formula`** — polarization identity: ⟨b−a, c−a⟩ = (d(a,b)² + d(a,c)² − d(b,c)²)/2.
6. **`gram3_entry_00` through `gram3_entry_22`** — all 9 Gram matrix entries expressed in terms of squared distances.
7. **`gram3_eq_gramFromDistances`** — the Gram matrix equals the explicit `gramFromDistances` construction from the six pairwise squared distances.
8. **`orient_sq_determined_by_pairwise_distances`** — orient² = det(gramFromDistances(d_ab, d_ac, d_ad, d_bc, d_bd, d_cd)).
9. **`same_distances_same_orient_sq_opposite_orient_standardTet`** — the standard tetrahedron and its z-reflection have equal distance matrices and equal orient², but opposite orient.
10. **`distance_data_determines_magnitude_not_chirality_sign`** — same distances and orient², but orient values differ (≠).
11. **`no_distance_only_orient_recoverer_standardTet`** — no function R from distance matrices to ℝ can recover the true oriented volume for both a tetrahedron and its mirror image.
12. **`no_distance_only_sign_recoverer_standardTet`** — no Bool-valued classifier can separate a tetrahedron from its mirror image using only distances.

The file imports and reuses definitions from `ChiralityNoGo.lean` and `OrthogonalChirality.lean`.