/-
# Orthogonal Transformations, Distances, and Chirality

Extends ChiralityNoGo.lean to the general orthogonal case:
- Pairwise distances are invariant under orthogonal transformations.
- Oriented volume transforms by the determinant of the transformation.
- Improper rotations (det = -1) flip chirality; proper rotations (det = 1) preserve it.
-/
import Mathlib
import RequestProject.ChiralityNoGo

open Finset BigOperators Matrix

-- ============================================================
-- SECTION 1. Definitions
-- ============================================================

/-- Matrix-vector multiplication for 3×3 real matrices and R³ vectors. -/
def matVec (A : Matrix (Fin 3) (Fin 3) ℝ) (p : Point3) : Point3 :=
  fun i => ∑ j : Fin 3, A i j * p j

/-- Transform a point by a 3×3 matrix. -/
def transformPoint (A : Matrix (Fin 3) (Fin 3) ℝ) (p : Point3) : Point3 :=
  matVec A p

/-- Transform all four points in a configuration. -/
def transformConfig4 (A : Matrix (Fin 3) (Fin 3) ℝ) (P : Config4) : Config4 :=
  fun k => transformPoint A (P k)

/-- A 3×3 real matrix is orthogonal iff AᵀA = I. -/
def IsOrthogonal3 (A : Matrix (Fin 3) (Fin 3) ℝ) : Prop :=
  Aᵀ * A = 1

/-- The edge matrix of four points: columns are b-a, c-a, d-a. -/
def edgeMatrix (a b c d : Point3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun r j => colVec a b c d j r

-- ============================================================
-- SECTION 2. Key algebraic lemmas
-- ============================================================

/-- transformPoint is the same as Matrix.mulVec. -/
theorem transformPoint_eq_mulVec (A : Matrix (Fin 3) (Fin 3) ℝ) (p : Point3) :
    transformPoint A p = A.mulVec p := by
  ext i; simp [transformPoint, matVec, Matrix.mulVec, dotProduct]

/-- transformPoint distributes: A(p - q) componentwise. -/
theorem transformPoint_sub (A : Matrix (Fin 3) (Fin 3) ℝ) (p q : Point3) (i : Fin 3) :
    transformPoint A p i - transformPoint A q i = ∑ j : Fin 3, A i j * (p j - q j) := by
  simp only [transformPoint, matVec]
  rw [← Finset.sum_sub_distrib]
  congr 1; funext j; ring

/-
Orthogonality gives AᵀA = I entrywise.
-/
theorem orthogonal_entry (A : Matrix (Fin 3) (Fin 3) ℝ) (hA : IsOrthogonal3 A)
    (j k : Fin 3) :
    ∑ i : Fin 3, A i j * A i k = if j = k then 1 else 0 := by
  replace hA := congr_fun ( congr_fun hA j ) k; fin_cases j <;> fin_cases k <;> simp_all +decide [ Matrix.mul_apply ] ;

/-
distSq expressed via orthogonality.
-/
theorem distSq_transformPoint_of_orthogonal
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (p q : Point3) :
    distSq (transformPoint A p) (transformPoint A q) = distSq p q := by
  unfold distSq;
  simp +decide only [transformPoint_sub, mul_sub, pow_two];
  have := orthogonal_entry A hA; simp_all +decide [Fin.sum_univ_three, mul_sub, mul_comm];
  grind +ring

-- ============================================================
-- SECTION 3. Distance matrix invariance
-- ============================================================

/-- The full pairwise distance matrix is invariant under orthogonal transformations. -/
theorem distMatrix4_transformConfig4_of_orthogonal
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (P : Config4) :
    distMatrix4 (transformConfig4 A P) = distMatrix4 P := by
  ext i j
  exact distSq_transformPoint_of_orthogonal A hA (P i) (P j)

-- ============================================================
-- SECTION 4. Orientation transforms by the determinant
-- ============================================================

/-- orient is det of the edge matrix. -/
theorem orient_eq_det_edgeMatrix (a b c d : Point3) :
    orient a b c d = Matrix.det (edgeMatrix a b c d) := by
  rfl

/-
The edge matrix of transformed points equals A times the original edge matrix.
-/
theorem edgeMatrix_transform (A : Matrix (Fin 3) (Fin 3) ℝ) (a b c d : Point3) :
    edgeMatrix (transformPoint A a) (transformPoint A b)
               (transformPoint A c) (transformPoint A d) =
    A * edgeMatrix a b c d := by
  unfold edgeMatrix; simp +decide [ transformPoint ] ;
  unfold colVec matVec; ext r j; fin_cases r <;> fin_cases j <;> simp +decide [ Matrix.mul_apply, Fin.sum_univ_three ] <;> ring;

/-- Oriented volume transforms by det(A). -/
theorem orient_transformPoint
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (a b c d : Point3) :
    orient (transformPoint A a) (transformPoint A b)
           (transformPoint A c) (transformPoint A d) =
    Matrix.det A * orient a b c d := by
  rw [orient_eq_det_edgeMatrix, orient_eq_det_edgeMatrix]
  rw [edgeMatrix_transform]
  rw [Matrix.det_mul]

-- ============================================================
-- SECTION 5. Consequences for chirality
-- ============================================================

/-- If det(A) = -1, orientation is negated. -/
theorem orientation_flips_of_det_neg_one
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hdet : Matrix.det A = -1)
    (a b c d : Point3) :
    orient (transformPoint A a) (transformPoint A b)
           (transformPoint A c) (transformPoint A d) =
    -orient a b c d := by
  rw [orient_transformPoint, hdet]; ring

/-- If det(A) = 1, orientation is preserved. -/
theorem orientation_preserved_of_det_one
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hdet : Matrix.det A = 1)
    (a b c d : Point3) :
    orient (transformPoint A a) (transformPoint A b)
           (transformPoint A c) (transformPoint A d) =
    orient a b c d := by
  rw [orient_transformPoint, hdet]; ring

-- ============================================================
-- SECTION 6. General impossibility and obstruction theorems
-- ============================================================

/-- Any function of the distance matrix is invariant under orthogonal transformations. -/
theorem distance_only_function_orthogonal_invariant
    {α : Type}
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (F : Matrix (Fin 4) (Fin 4) ℝ → α)
    (P : Config4) :
    F (distMatrix4 (transformConfig4 A P)) = F (distMatrix4 P) := by
  rw [distMatrix4_transformConfig4_of_orthogonal A hA P]

/-- Improper orthogonal transformations preserve distances but flip chirality. -/
theorem improper_orthogonal_distance_chirality_obstruction
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (hdet : Matrix.det A = -1)
    (P : Config4) :
    distMatrix4 (transformConfig4 A P) = distMatrix4 P ∧
    orient4 (transformConfig4 A P) = -orient4 P := by
  constructor
  · exact distMatrix4_transformConfig4_of_orthogonal A hA P
  · unfold orient4 transformConfig4
    exact orientation_flips_of_det_neg_one A hdet (P 0) (P 1) (P 2) (P 3)

/-- Proper orthogonal transformations preserve both distances and chirality. -/
theorem proper_orthogonal_preserves_all
    (A : Matrix (Fin 3) (Fin 3) ℝ)
    (hA : IsOrthogonal3 A)
    (hdet : Matrix.det A = 1)
    (P : Config4) :
    distMatrix4 (transformConfig4 A P) = distMatrix4 P ∧
    orient4 (transformConfig4 A P) = orient4 P := by
  constructor
  · exact distMatrix4_transformConfig4_of_orthogonal A hA P
  · unfold orient4 transformConfig4
    exact orientation_preserved_of_det_one A hdet (P 0) (P 1) (P 2) (P 3)