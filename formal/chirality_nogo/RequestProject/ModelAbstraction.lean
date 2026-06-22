/-
# Model-Level Abstraction: Distance-Factoring Models Cannot Represent Chirality

Any model whose geometric input factors through the pairwise distance matrix
is reflection-invariant, and therefore cannot represent any target function
that changes sign under reflection.
-/
import Mathlib
import RequestProject.ChiralityNoGo

open Finset BigOperators Matrix

-- ============================================================
-- SECTION 1. FactorsThroughDistMatrix
-- ============================================================

/-- A model `M : Config4 → α` factors through the distance matrix if there exists
    `F : Matrix (Fin 4) (Fin 4) ℝ → α` such that `M P = F (distMatrix4 P)` for all `P`. -/
def FactorsThroughDistMatrix
    {α : Type}
    (M : Config4 → α) : Prop :=
  ∃ F : Matrix (Fin 4) (Fin 4) ℝ → α,
    ∀ P : Config4, M P = F (distMatrix4 P)

-- ============================================================
-- SECTION 2. Reflection invariance of distance-factoring models
-- ============================================================

/-- Any model that factors through the distance matrix is reflection-invariant. -/
theorem factorsThroughDistMatrix_reflection_invariant
    {α : Type}
    (M : Config4 → α)
    (hM : FactorsThroughDistMatrix M)
    (P : Config4) :
    M (reflectZConfig P) = M P := by
  obtain ⟨F, hF⟩ := hM
  rw [hF, hF, distMatrix4_reflectZConfig]

-- ============================================================
-- SECTION 3. orientTarget impossibility
-- ============================================================

/-- The orientation target function. -/
def orientTarget (P : Config4) : ℝ := orient4 P

/-- `orientTarget` does not factor through the distance matrix. -/
theorem orientTarget_not_factorsThroughDistMatrix :
    ¬ FactorsThroughDistMatrix orientTarget := by
  intro hM
  have hinv := factorsThroughDistMatrix_reflection_invariant orientTarget hM standardTet
  unfold orientTarget at hinv
  rw [orient4_reflectZConfig] at hinv
  have h1 := orient4_standardTet
  linarith

-- ============================================================
-- SECTION 4. General reflection-odd targets
-- ============================================================

/-- Any target that flips sign under reflection at some nondegenerate configuration
    cannot factor through the distance matrix. -/
theorem reflection_odd_target_not_distance_factorable
    (T : Config4 → ℝ)
    (P : Config4)
    (hflip : T (reflectZConfig P) = -T P)
    (hnonzero : T P ≠ 0) :
    ¬ FactorsThroughDistMatrix T := by
  intro hM
  have hinv := factorsThroughDistMatrix_reflection_invariant T hM P
  have : T P = -T P := by linarith
  have : 2 * T P = 0 := by linarith
  have : T P = 0 := by linarith
  exact hnonzero this

/-- `orient4` does not factor through the distance matrix (via the general theorem). -/
theorem orient4_not_distance_factorable_general :
    ¬ FactorsThroughDistMatrix orient4 :=
  reflection_odd_target_not_distance_factorable
    orient4
    standardTet
    (by rw [orient4_reflectZConfig, orient4_standardTet])
    (by rw [orient4_standardTet]; exact one_ne_zero)


-- ============================================================
-- SECTION 5. Bool chirality sign version
-- ============================================================

/-- Boolean chirality sign: true iff orient4 P > 0. -/
noncomputable def chiralitySignBool (P : Config4) : Bool :=
  if 0 < orient4 P then true else false

theorem chiralitySignBool_standardTet :
    chiralitySignBool standardTet = true := by
  unfold chiralitySignBool
  rw [orient4_standardTet]
  norm_num

theorem chiralitySignBool_reflect_standardTet :
    chiralitySignBool (reflectZConfig standardTet) = false := by
  unfold chiralitySignBool
  rw [orient4_reflectZ_standardTet]
  norm_num

/-- The Boolean chirality sign does not factor through the distance matrix. -/
theorem chiralitySignBool_not_distance_factorable :
    ¬ FactorsThroughDistMatrix chiralitySignBool := by
  intro ⟨F, hF⟩
  have h1 := hF standardTet
  have h2 := hF (reflectZConfig standardTet)
  rw [distMatrix4_reflectZConfig] at h2
  rw [chiralitySignBool_standardTet] at h1
  rw [chiralitySignBool_reflect_standardTet] at h2
  simp_all
