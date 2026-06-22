/-
# Distances Determine Squared Oriented Volume, But Not the Sign

This file proves that pairwise squared distances among four points in ℝ³
determine the square of the oriented volume (via the Gram matrix determinant),
but cannot determine its sign.
-/
import Mathlib
import RequestProject.ChiralityNoGo
import RequestProject.OrthogonalChirality

open Finset BigOperators Matrix

-- ============================================================
-- SECTION 1. Setup: vector subtraction, dot product, norm
-- ============================================================

/-- Vector subtraction of two points. -/
def vsub (p q : Point3) : Point3 :=
  fun i => p i - q i

/-- Dot product on ℝ³. -/
def dot3 (x y : Point3) : ℝ :=
  ∑ i : Fin 3, x i * y i

/-- Squared norm. -/
def normSq3 (x : Point3) : ℝ :=
  dot3 x x

theorem distSq_eq_normSq_vsub (p q : Point3) :
    distSq p q = normSq3 (vsub p q) := by
  simp [distSq, normSq3, dot3, vsub, sq, sub_mul]

-- ============================================================
-- SECTION 2. Edge vectors and Gram matrix
-- ============================================================

/-- Edge vector u = b - a. -/
def edgeU (a b _c _d : Point3) : Point3 := vsub b a
/-- Edge vector v = c - a. -/
def edgeV (a _b c _d : Point3) : Point3 := vsub c a
/-- Edge vector w = d - a. -/
def edgeW (a _b _c d : Point3) : Point3 := vsub d a

/-- Edge matrix whose columns are u, v, w. -/
def edgeMatrixFromPoints (a b c d : Point3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j =>
    if j = 0 then edgeU a b c d i
    else if j = 1 then edgeV a b c d i
    else edgeW a b c d i

/-- Gram matrix: G_{ij} = dot(column_i, column_j) of the edge matrix. -/
def gram3 (a b c d : Point3) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j =>
    dot3
      (fun k => edgeMatrixFromPoints a b c d k i)
      (fun k => edgeMatrixFromPoints a b c d k j)

theorem gram3_eq_transpose_mul_edgeMatrix
    (a b c d : Point3) :
    gram3 a b c d =
      (edgeMatrixFromPoints a b c d)ᵀ * edgeMatrixFromPoints a b c d := by
  ext i j
  simp [gram3, dot3, Matrix.mul_apply, Matrix.transpose_apply]

-- ============================================================
-- SECTION 3. Determinant-square theorem
-- ============================================================

theorem orient_eq_det_edgeMatrixFromPoints
    (a b c d : Point3) :
    orient a b c d = Matrix.det (edgeMatrixFromPoints a b c d) := by
  unfold orient edgeMatrixFromPoints edgeU edgeV edgeW vsub colVec
  congr 1
  ext i j
  fin_cases j <;> simp

theorem det_gram3_eq_orient_sq
    (a b c d : Point3) :
    Matrix.det (gram3 a b c d) = orient a b c d ^ 2 := by
  rw [gram3_eq_transpose_mul_edgeMatrix, Matrix.det_mul, Matrix.det_transpose,
      orient_eq_det_edgeMatrixFromPoints]
  ring

-- ============================================================
-- SECTION 4. Distances determine Gram entries (polarization)
-- ============================================================

theorem dot_vsub_vsub_eq_distSq_formula
    (a b c : Point3) :
    dot3 (vsub b a) (vsub c a) =
      (distSq a b + distSq a c - distSq b c) / 2 := by
  simp [dot3, distSq, vsub, Fin.sum_univ_three]
  ring

theorem gram3_entry_00 (a b c d : Point3) :
    gram3 a b c d 0 0 = distSq a b := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeU, vsub, distSq, Fin.sum_univ_three, sq]
  ring

theorem gram3_entry_11 (a b c d : Point3) :
    gram3 a b c d 1 1 = distSq a c := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeV, vsub, distSq, Fin.sum_univ_three, sq]
  ring

theorem gram3_entry_22 (a b c d : Point3) :
    gram3 a b c d 2 2 = distSq a d := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeW, vsub, distSq, Fin.sum_univ_three, sq]
  ring

theorem gram3_entry_01 (a b c d : Point3) :
    gram3 a b c d 0 1 =
      (distSq a b + distSq a c - distSq b c) / 2 := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeU, edgeV, vsub, distSq, Fin.sum_univ_three, sq]
  ring

theorem gram3_entry_02 (a b c d : Point3) :
    gram3 a b c d 0 2 =
      (distSq a b + distSq a d - distSq b d) / 2 := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeU, edgeW, vsub, distSq, Fin.sum_univ_three, sq]
  ring

theorem gram3_entry_12 (a b c d : Point3) :
    gram3 a b c d 1 2 =
      (distSq a c + distSq a d - distSq c d) / 2 := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeV, edgeW, vsub, distSq, Fin.sum_univ_three, sq]
  ring

theorem gram3_entry_10 (a b c d : Point3) :
    gram3 a b c d 1 0 =
      (distSq a b + distSq a c - distSq b c) / 2 := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeU, edgeV, vsub, distSq, Fin.sum_univ_three, sq]
  ring

theorem gram3_entry_20 (a b c d : Point3) :
    gram3 a b c d 2 0 =
      (distSq a b + distSq a d - distSq b d) / 2 := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeU, edgeW, vsub, distSq, Fin.sum_univ_three, sq]
  ring

theorem gram3_entry_21 (a b c d : Point3) :
    gram3 a b c d 2 1 =
      (distSq a c + distSq a d - distSq c d) / 2 := by
  simp +decide [gram3, dot3, edgeMatrixFromPoints, edgeV, edgeW, vsub, distSq, Fin.sum_univ_three, sq]
  ring

-- ============================================================
-- SECTION 5. Explicit Gram-from-distance matrix
-- ============================================================

/-- Gram matrix constructed directly from six pairwise squared distances. -/
noncomputable def gramFromDistances
    (dab dac dad dbc dbd dcd : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  fun i j =>
    if i = 0 ∧ j = 0 then dab
    else if i = 1 ∧ j = 1 then dac
    else if i = 2 ∧ j = 2 then dad
    else if (i = 0 ∧ j = 1) ∨ (i = 1 ∧ j = 0) then (dab + dac - dbc) / 2
    else if (i = 0 ∧ j = 2) ∨ (i = 2 ∧ j = 0) then (dab + dad - dbd) / 2
    else (dac + dad - dcd) / 2

theorem gram3_eq_gramFromDistances
    (a b c d : Point3) :
    gram3 a b c d =
      gramFromDistances
        (distSq a b)
        (distSq a c)
        (distSq a d)
        (distSq b c)
        (distSq b d)
        (distSq c d) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp +decide [gramFromDistances,
      gram3_entry_00, gram3_entry_01, gram3_entry_02,
      gram3_entry_10, gram3_entry_11, gram3_entry_12,
      gram3_entry_20, gram3_entry_21, gram3_entry_22]

theorem orient_sq_determined_by_pairwise_distances
    (a b c d : Point3) :
    orient a b c d ^ 2 =
      Matrix.det
        (gramFromDistances
          (distSq a b)
          (distSq a c)
          (distSq a d)
          (distSq b c)
          (distSq b d)
          (distSq c d)) := by
  rw [← det_gram3_eq_orient_sq, gram3_eq_gramFromDistances]

-- ============================================================
-- SECTION 6. Sign is not determined
-- ============================================================

theorem same_distances_same_orient_sq_opposite_orient_standardTet :
    distMatrix4 (reflectZConfig standardTet) = distMatrix4 standardTet ∧
    orient4 (reflectZConfig standardTet) ^ 2 = orient4 standardTet ^ 2 ∧
    orient4 (reflectZConfig standardTet) = -orient4 standardTet := by
  refine ⟨distMatrix4_standardTet_reflect, ?_, orient4_reflectZConfig standardTet⟩
  rw [orient4_reflectZ_standardTet, orient4_standardTet]
  norm_num

theorem distance_data_determines_magnitude_not_chirality_sign :
    distMatrix4 (reflectZConfig standardTet) = distMatrix4 standardTet ∧
    orient4 (reflectZConfig standardTet) ^ 2 = orient4 standardTet ^ 2 ∧
    orient4 (reflectZConfig standardTet) ≠ orient4 standardTet := by
  refine ⟨distMatrix4_standardTet_reflect, ?_, ?_⟩
  · rw [orient4_reflectZ_standardTet, orient4_standardTet]; norm_num
  · rw [orient4_reflectZ_standardTet, orient4_standardTet]; norm_num

-- ============================================================
-- SECTION 7. Classifier impossibility
-- ============================================================

/-- A sign recoverer is any function from the distance matrix to ℝ. -/
abbrev SignRecoverer := Matrix (Fin 4) (Fin 4) ℝ → ℝ

theorem no_distance_only_orient_recoverer_standardTet
    (R : Matrix (Fin 4) (Fin 4) ℝ → ℝ) :
    ¬ (
      R (distMatrix4 standardTet) = orient4 standardTet ∧
      R (distMatrix4 (reflectZConfig standardTet)) =
        orient4 (reflectZConfig standardTet)
    ) := by
  rw [distMatrix4_standardTet_reflect]
  intro ⟨h1, h2⟩
  rw [h1] at h2
  rw [orient4_standardTet, orient4_reflectZ_standardTet] at h2
  norm_num at h2

theorem no_distance_only_sign_recoverer_standardTet
    (R : Matrix (Fin 4) (Fin 4) ℝ → Bool) :
    ¬ (
      R (distMatrix4 standardTet) = true ∧
      R (distMatrix4 (reflectZConfig standardTet)) = false
    ) := by
  rw [distMatrix4_standardTet_reflect]
  intro ⟨h1, h2⟩
  rw [h1] at h2
  exact absurd h2 (by decide)
