# Harmonic Aristotle — ChiralFold Lean 4 Proof Generalization Prompt

**Purpose:** Paste this entire document into [Harmonic Aristotle](https://aristotle.harmonic.fun/) (or your Aristotle workspace) together with the Lean 4 project at `formal/chirality_nogo/` from the ChiralFold repository (`cursor/aristotle-formal-proofs-9901` branch or `formal/chirality_nogo/` after merge).

**Goal:** Generalize and strengthen the existing no-go theorems so they address reviewer concerns that the current proofs only apply to naive “distance-matrix-only” models and not to SE(3)-equivariant neural networks (e.g. AlphaFold 3).

---

## Context for Aristotle

ChiralFold is a protein stereochemistry toolkit. We formalized in Lean 4 that **pairwise distances among four points in ℝ³ determine the squared oriented volume (Gram determinant) but not its sign**, and that **any model factoring through the distance matrix cannot represent signed chirality**.

The manuscript uses this to motivate why diffusion predictors need explicit signed-volume (or equivalent) chirality checks—not to claim a complete impossibility theorem about all neural networks.

### Existing Lean 4 project structure

```
formal/chirality_nogo/
  lakefile.toml
  lean-toolchain          # Lean 4.28.0
  RequestProject/
    Main.lean
    ChiralityNoGo.lean       # Core: reflection preserves distMatrix4, negates orient
    OrthogonalChirality.lean # Orthogonal maps preserve distances; det flips orient
    ModelAbstraction.lean    # FactorsThroughDistMatrix → reflection invariant
    DistancesOrientVolume.lean # det(Gram) = orient²; distances fix |orient| not sign
```

### Key definitions already formalized

- `Point3`, `Config4` — four points in ℝ³  
- `distSq`, `distMatrix4` — pairwise squared distances  
- `reflectZConfig` — reflection through xy-plane (improper isometry, det = −1)  
- `orient` / `orient4` — signed tetrahedron volume (Cα chirality surrogate)  
- `FactorsThroughDistMatrix M` — ∃ F, M P = F (distMatrix4 P)  
- `gram3`, `edgeMatrixFromPoints` — Gram matrix link to oriented volume  

### Key theorems already proved (machine-checked)

| Theorem | File | Statement (informal) |
|---------|------|----------------------|
| `distMatrix4_reflectZConfig` | ChiralityNoGo | Reflection preserves all pairwise squared distances |
| `orient_reflectZConfig` | ChiralityNoGo | Reflection negates oriented volume |
| `same_distances_opposite_orientations_standardTet` | ChiralityNoGo | Equal distance matrices, opposite orientations exist |
| `factorsThroughDistMatrix_reflection_invariant` | ModelAbstraction | Distance-factoring models are reflection-invariant |
| `orient4_not_distance_factorable_general` | ModelAbstraction | orient4 cannot factor through distMatrix4 |
| `chiralitySignBool_not_distance_factorable` | ModelAbstraction | Boolean chirality sign cannot factor through distances |
| `det_gram3_eq_orient_sq` | DistancesOrientVolume | det(Gram) = orient² |
| `distance_data_determines_magnitude_not_chirality_sign` | DistancesOrientVolume | Distances determine \|orient\|² but not sign |

Build: `cd formal/chirality_nogo && lake build`

---

## Task for Aristotle

Please **generalize, tighten, and extend** this formalization. Work incrementally; preserve all existing theorems unless you find a bug. Target Mathlib 4.28 compatibility.

### Priority 1 — Strengthen the abstraction layer

1. **Generalize beyond `reflectZConfig`** to any **improper orthogonal** transformation (orthogonal with det = −1), not just z-reflection. The file `OrthogonalChirality.lean` begins this; complete it if gaps remain.

2. Prove a clean **theorem statement** of the form:

   > If `M : Config4 → α` is invariant under all improper isometries of ℝ³ that fix the distance matrix, and chirality sign flips under such an isometry, then `M` cannot factor through `distMatrix4`.

3. Add a corollary for **Cα chirality in proteins**:

   > The signed volume of the tetrahedron (N, Cα, C, Cβ) cannot be a function of the six pairwise distances among those four atoms alone.

   Use the existing `orient4` API or define `proteinCaOrient` with an explicit connection to the four-point configuration.

### Priority 2 — Bridge to “equivariant models” (careful scope)

We do **not** want an overclaim about all neural networks. We want precise intermediate results:

4. Define a class `DistanceOnlyRepresentation` (or similar): maps `Config4 → ℝⁿ` that factors through `distMatrix4` **and** is continuous (or smooth) if helpful.

5. Prove: **no continuous classifier** `Config4 → Bool` that factors through distances can equal `chiralitySignBool` on all configurations.

6. Optional: formalize that **SE(3)-equivariant maps whose output is invariant under the full distance matrix** (i.e. blind to orientation of the volume form) cannot distinguish enantiomers. State hypotheses explicitly—do not claim this covers all GNNs with angular features.

7. Add a short `README.md` in `formal/chirality_nogo/` with:
   - What is proved (distance-only obstruction)
   - What is **not** proved (full AF3 architecture)
   - Suggested one-sentence manuscript wording

### Priority 3 — Paper-facing exports

8. Produce a **LaTeX snippet** (plain text) suitable for Bioinformatics Methods §2.8:

   - One paragraph describing the formal result
   - One table mapping theorem names → biological consequence
   - Explicit limitation sentence about SE(3)-equivariant predictors

9. Add `RequestProject/Examples.lean` with:
   - `standardTet` and its mirror configuration
   - `norm_num` proofs that distance matrices match but orientations differ

### Priority 4 — Quality bar

- All new proofs must compile with `lake build` with no `sorry`.
- Prefer reusable lemmas over one-off `simp` chains.
- If a conjecture is too hard (e.g. full EGNN obstruction), state it as `theorem ... := by sorry` **only** in a separate `Conjectures.lean` file clearly marked not for publication.
- Keep `maxHeartbeats` reasonable; split heavy proofs across lemmas.

---

## Reviewer concern to address (quote for Aristotle)

> “The Lean proofs only show that chirality sign cannot be recovered from pairwise distances. AlphaFold 3 is not a distance-matrix model—it uses SE(3)-equivariant networks with angular features. The formal section overclaims relevance to AF3.”

**Desired response in math:** A hierarchy of results:

1. **Proved:** Distance data → \|orient\|² only; sign requires non-distance information (signed volume, torsion angles, or chirality-aware features).  
2. **Proved:** Any model whose *entire* input representation is the distance matrix cannot classify chirality.  
3. **Not claimed:** Impossibility for all equivariant predictors.  
4. **Suggested manuscript fix:** “Motivates explicit signed-volume auditing” rather than “explains AF3 failure.”

---

## Suggested manuscript wording (for Aristotle to refine)

> Machine-checked Lean 4 proofs establish that pairwise distances among four backbone atoms determine the magnitude but not the sign of the Cα signed tetrahedron volume, and that any representation factoring solely through the pairwise distance matrix is blind to enantiomeric configuration. These results characterize a structural limitation of distance-only geometric representations; they do not constitute a complete model of SE(3)-equivariant predictors such as AlphaFold 3, which incorporate angular and higher-order features. They motivate ChiralFold’s use of signed volume as a minimal non-distance invariant for deposition auditing.

---

## Files to attach when prompting Aristotle

1. Entire `formal/chirality_nogo/` directory (or tarball from branch `cursor/aristotle-formal-proofs-9901`)
2. This prompt (`ARISTOTLE_PROMPT.md`)
3. Optional: `paper/chiralfold_paper.tex` Methods §2.8 and Discussion paragraphs on formal verification

---

## Success criteria

- [ ] `lake build` passes with zero `sorry` in production files  
- [ ] At least one new theorem explicitly about **improper orthogonal** maps (not just z-reflection)  
- [ ] `README.md` in formal/ with proved / not-proved scope  
- [ ] LaTeX snippet ready to paste into submission  
- [ ] No regression of existing theorem names used in the paper supplementary table

---

## After Aristotle completes

1. Run `lake build` locally and in CI  
2. Update `paper/submission/bioinformatics/supplementary_material.tex` Table S2 if theorem names change  
3. Replace Discussion language per README scope section  
4. Commit to `formal/chirality_nogo/` and tag release
