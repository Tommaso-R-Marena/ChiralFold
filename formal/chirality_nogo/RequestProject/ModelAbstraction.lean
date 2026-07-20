/-
# Model-Level Abstraction: Distance-Factoring Models Cannot Represent Chirality

This file isolates the exact scope of the obstruction: it applies to models whose
entire geometric representation factors through the four-point distance matrix.
No claim is made about arbitrary SE(3)-equivariant architectures.
-/
import Mathlib
import RequestProject.ChiralityNoGo
import RequestProject.OrthogonalChirality

open Finset BigOperators Matrix

/-- A model factors through the distance matrix when all of its dependence on a
configuration passes through the matrix of pairwise squared distances. -/
def FactorsThroughDistMatrix
    {α : Type}
    (M : Config4 → α) : Prop :=
  ∃ F : Matrix (Fin 4) (Fin 4) ℝ → α,
    ∀ P : Config4, M P = F (distMatrix4 P)

/-
Any distance-factoring model takes equal values on configurations having the
same distance matrix. This is the reusable abstraction behind all no-go results.
-/
theorem factorsThroughDistMatrix_of_distMatrix_eq
    {α : Type}
    (M : Config4 → α)
    (hM : FactorsThroughDistMatrix M)
    {P Q : Config4}
    (hPQ : distMatrix4 P = distMatrix4 Q) :
    M P = M Q := by
  obtain ⟨ F, hF ⟩ := hM; aesop;

/-
Any model that factors through the distance matrix is reflection-invariant.
-/
theorem factorsThroughDistMatrix_reflection_invariant
    {α : Type}
    (M : Config4 → α)
    (hM : FactorsThroughDistMatrix M)
    (P : Config4) :
    M (reflectZConfig P) = M P := by
  obtain ⟨ F, hF ⟩ := hM;
  rw [ hF, hF, distMatrix4_reflectZConfig ]

/-
Distance-factoring models are invariant under every orthogonal coordinate
change, including every improper one.
-/
theorem factorsThroughDistMatrix_orthogonal_invariant
    {α : Type}
    (M : Config4 → α)
    (hM : FactorsThroughDistMatrix M)
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (P : Config4) :
    M (transformConfig4 A P) = M P := by
  obtain ⟨ F, hF ⟩ := hM;
  rw [ hF, hF, distMatrix4_transformConfig4_of_orthogonal A hA ]

/-
A clean witness theorem: if a distance-preserving transformation changes a
target at one configuration, that target cannot factor through distances.
-/
theorem distance_preserving_target_change_not_factorable
    {α : Type}
    (T : Config4 → α)
    (g : Config4 → Config4)
    (P : Config4)
    (hdist : distMatrix4 (g P) = distMatrix4 P)
    (hchange : T (g P) ≠ T P) :
    ¬ FactorsThroughDistMatrix T := by
  exact fun ⟨ F, hF ⟩ => hchange <| by simp +decide [ hF, hdist ] ;

/-
Specialization to an improper orthogonal map: if the chirality target changes
under such a map at one witness, it is not distance-factorable. The assumptions
state precisely the scope needed; no neural-network architecture is assumed.
-/
theorem improper_orthogonal_target_not_distance_factorable
    {α : Type}
    (T : Config4 → α)
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (_hdet : Matrix.det A = -1)
    (P : Config4)
    (hchange : T (transformConfig4 A P) ≠ T P) :
    ¬ FactorsThroughDistMatrix T := by
  exact fun h => hchange <| factorsThroughDistMatrix_orthogonal_invariant T h A hA P

/-
At every nondegenerate witness, any improper orthogonal coordinate change
exhibits directly that signed orientation is not distance-factorable.
-/
theorem orient4_not_distance_factorable_via_improper_orthogonal
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (hdet : Matrix.det A = -1)
    (P : Config4)
    (hnonzero : orient4 P ≠ 0) :
    ¬ FactorsThroughDistMatrix orient4 := by
  apply improper_orthogonal_target_not_distance_factorable orient4 A hA hdet P
  intro heq
  have hflip := (improper_orthogonal_distance_chirality_obstruction A hA hdet P).2
  apply hnonzero
  linarith

/-- The orientation target function. -/
def orientTarget (P : Config4) : ℝ := orient4 P

/-
`orientTarget` does not factor through the distance matrix.
-/
theorem orientTarget_not_factorsThroughDistMatrix :
    ¬ FactorsThroughDistMatrix orientTarget := by
  intro h
  have hinv := factorsThroughDistMatrix_reflection_invariant orientTarget h standardTet
  unfold orientTarget at hinv
  rw [orient4_reflectZ_standardTet, orient4_standardTet] at hinv
  norm_num at hinv

/-
Any real-valued target that flips sign under reflection at a nondegenerate
configuration cannot factor through the distance matrix.
-/
theorem reflection_odd_target_not_distance_factorable
    (T : Config4 → ℝ)
    (P : Config4)
    (hflip : T (reflectZConfig P) = -T P)
    (hnonzero : T P ≠ 0) :
    ¬ FactorsThroughDistMatrix T := by
  exact fun h => hnonzero <| by linarith [ factorsThroughDistMatrix_reflection_invariant T h P ] ;

/-
`orient4` does not factor through the distance matrix.
-/
theorem orient4_not_distance_factorable_general :
    ¬ FactorsThroughDistMatrix orient4 := by
  convert orientTarget_not_factorsThroughDistMatrix using 1

/-- Boolean chirality sign: true iff `orient4 P > 0`. -/
noncomputable def chiralitySignBool (P : Config4) : Bool :=
  if 0 < orient4 P then true else false

theorem chiralitySignBool_standardTet :
    chiralitySignBool standardTet = true := by
  unfold chiralitySignBool; norm_num [ orient4_standardTet ] ;

theorem chiralitySignBool_reflect_standardTet :
    chiralitySignBool (reflectZConfig standardTet) = false := by
  unfold chiralitySignBool;
  rw [ orient4_reflectZ_standardTet ] ; norm_num

/-
The Boolean chirality sign does not factor through the distance matrix.
-/
theorem chiralitySignBool_not_distance_factorable :
    ¬ FactorsThroughDistMatrix chiralitySignBool := by
  grind +suggestions

/-- A distance-only real-vector representation, with continuity recorded
separately from the algebraic factorization property. Continuity does not evade
the obstruction (indeed, the no-go theorem needs only factorization). -/
structure DistanceOnlyRepresentation (n : ℕ) where
  toFun : Config4 → (Fin n → ℝ)
  factors : FactorsThroughDistMatrix toFun
  continuous_toFun : Continuous toFun

instance (n : ℕ) : CoeFun (DistanceOnlyRepresentation n)
    (fun _ => Config4 → (Fin n → ℝ)) := ⟨DistanceOnlyRepresentation.toFun⟩

/-
No classifier—continuous or otherwise—that factors through distances can
agree everywhere with Boolean chirality sign. Here continuity of the real-vector
representation is explicit, while the final classifier may be arbitrary.
-/
theorem no_continuous_distance_representation_classifier
    (n : ℕ)
    (R : DistanceOnlyRepresentation n)
    (C : (Fin n → ℝ) → Bool) :
    ¬ ∀ P : Config4, C (R P) = chiralitySignBool P := by
  intro h;
  convert chiralitySignBool_not_distance_factorable _;
  exact ⟨ fun D => C ( R.factors.choose D ), fun P => h P ▸ by rw [ R.factors.choose_spec P ] ⟩

/-
Direct formulation requested for continuous classifiers. Since `Bool` has its
standard discrete topology, continuity is an additional hypothesis; distance
factorization alone already gives the contradiction.
-/
theorem no_continuous_distance_factorable_chirality_classifier
    (C : Config4 → Bool)
    (_hcontinuous : Continuous C)
    (hfactor : FactorsThroughDistMatrix C) :
    ¬ ∀ P : Config4, C P = chiralitySignBool P := by
  obtain ⟨ F, hF ⟩ := hfactor;
  exact fun h => chiralitySignBool_not_distance_factorable ⟨ F, fun P => hF P ▸ h P ▸ rfl ⟩

/-
Explicitly scoped equivariant-model bridge. A model whose *output* is blind
to every orthogonal change of frame cannot distinguish an improper image from
the original. This does not assert that all SE(3)-equivariant networks satisfy
`hFullOrthogonalInvariant`.
-/
theorem full_orthogonal_invariant_model_cannot_distinguish_enantiomers
    {α : Type}
    (M : Config4 → α)
    (hFullOrthogonalInvariant :
      ∀ (A : Matrix (Fin 3) (Fin 3) ℝ), IsOrthogonal3 A →
        ∀ P, M (transformConfig4 A P) = M P)
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (_hdet : Matrix.det A = -1)
    (P : Config4) :
    M (transformConfig4 A P) = M P := by
  exact hFullOrthogonalInvariant A hA P

/-- Cα signed chirality for atoms ordered as `(N, Cα, C, Cβ)`. -/
def proteinCaOrient (n ca c cb : Point3) : ℝ := orient n ca c cb

/-- The configuration-level Cα target, with the atom order made explicit. -/
def proteinCaOrient4 (P : Config4) : ℝ :=
  proteinCaOrient (P 0) (P 1) (P 2) (P 3)

theorem proteinCaOrient4_eq_orient4 (P : Config4) :
    proteinCaOrient4 P = orient4 P := by
  rfl

/-
The signed volume of `(N, Cα, C, Cβ)` cannot be a function of the six
pairwise distances among those atoms alone.
-/
theorem proteinCaOrient_not_function_of_six_pairwise_distances :
    ¬ FactorsThroughDistMatrix proteinCaOrient4 := by
  convert orient4_not_distance_factorable_general using 1