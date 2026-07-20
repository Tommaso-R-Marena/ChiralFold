/-
# Pairwise Distances Do Not Determine Chirality

A formal no-go theorem: reflecting a tetrahedral configuration in R³ preserves
every pairwise squared distance but negates the signed oriented volume.
Therefore no function of the pairwise distance matrix alone can distinguish
a chiral center from its mirror image.
-/
import Mathlib

-- ============================================================
-- SECTION 1. Basic definitions
-- ============================================================

/-- A point in ℝ³. -/
abbrev Point3 := Fin 3 → ℝ

/-- A 4-point configuration in ℝ³. -/
abbrev Config4 := Fin 4 → Point3

/-- Squared Euclidean distance (purely algebraic). -/
def distSq (p q : Point3) : ℝ :=
  ∑ i : Fin 3, (p i - q i) ^ 2

/-- Pairwise squared distance matrix of a 4-point configuration. -/
def distMatrix4 (P : Config4) : Matrix (Fin 4) (Fin 4) ℝ :=
  fun i j => distSq (P i) (P j)

/-- Reflection through the xy-plane: negate the z-coordinate. -/
def reflectZPoint (p : Point3) : Point3 :=
  fun i => if i = (2 : Fin 3) then -p i else p i

/-- Reflect every point in a 4-point configuration. -/
def reflectZConfig (P : Config4) : Config4 :=
  fun k => reflectZPoint (P k)

-- ============================================================
-- SECTION 2. Oriented volume / chirality
-- ============================================================

/-- Column vector helper for the orientation matrix. -/
def colVec (a b c d : Point3) (j : Fin 3) : Point3 :=
  if j = 0 then fun r => b r - a r
  else if j = 1 then fun r => c r - a r
  else fun r => d r - a r

/-- Oriented volume determinant: det [b-a | c-a | d-a]. -/
def orient (a b c d : Point3) : ℝ :=
  Matrix.det (fun r j => colVec a b c d j r)

/-- Orientation of a 4-point configuration. -/
def orient4 (P : Config4) : ℝ :=
  orient (P 0) (P 1) (P 2) (P 3)

-- ============================================================
-- SECTION 3. Core algebraic lemmas
-- ============================================================

theorem reflectZPoint_involutive (p : Point3) :
    reflectZPoint (reflectZPoint p) = p := by
  funext i
  simp [reflectZPoint]
  split <;> ring

theorem distSq_reflectZPoint (p q : Point3) :
    distSq (reflectZPoint p) (reflectZPoint q) = distSq p q := by
  simp only [distSq, reflectZPoint]
  congr 1
  funext i
  fin_cases i
  · simp
  · simp
  · simp
    ring

theorem distMatrix4_reflectZConfig (P : Config4) :
    distMatrix4 (reflectZConfig P) = distMatrix4 P := by
  ext i j
  exact distSq_reflectZPoint (P i) (P j)

theorem orient_reflectZPoint (a b c d : Point3) :
    orient (reflectZPoint a) (reflectZPoint b) (reflectZPoint c) (reflectZPoint d) =
      -orient a b c d := by
  unfold orient colVec;
  simp +decide [ Matrix.det_fin_three, reflectZPoint ] ; ring;

theorem orient4_reflectZConfig (P : Config4) :
    orient4 (reflectZConfig P) = -orient4 P := by
  simp only [orient4, reflectZConfig]
  exact orient_reflectZPoint (P 0) (P 1) (P 2) (P 3)

-- ============================================================
-- SECTION 4. Standard tetrahedron witness
-- ============================================================

def e0 : Point3 := fun _ => 0
def e1 : Point3 := fun i => if i = 0 then 1 else 0
def e2 : Point3 := fun i => if i = 1 then 1 else 0
def e3 : Point3 := fun i => if i = 2 then 1 else 0

def standardTet : Config4 :=
  fun k =>
    if k = 0 then e0
    else if k = 1 then e1
    else if k = 2 then e2
    else e3

theorem orient4_standardTet : orient4 standardTet = 1 := by
  unfold orient4 standardTet;
  unfold orient e0 e1 e2 e3;
  -- The matrix is the 3x3 identity matrix, so its determinant is 1.
  simp [Matrix.det_fin_three];
  simp +decide [ colVec ]

theorem orient4_reflectZ_standardTet :
    orient4 (reflectZConfig standardTet) = -1 := by
  rw [orient4_reflectZConfig, orient4_standardTet]

theorem distMatrix4_standardTet_reflect :
    distMatrix4 (reflectZConfig standardTet) = distMatrix4 standardTet :=
  distMatrix4_reflectZConfig standardTet

-- ============================================================
-- SECTION 5. Distance-only classifier impossibility
-- ============================================================

theorem no_distance_only_classifier_separates_standardTet
    (C : Matrix (Fin 4) (Fin 4) ℝ → Bool) :
    ¬(C (distMatrix4 standardTet) = true ∧
      C (distMatrix4 (reflectZConfig standardTet)) = false) := by
  rw [distMatrix4_standardTet_reflect]
  intro ⟨h1, h2⟩
  simp_all

theorem no_distance_only_classifier_separates_standardTet_symmetric
    (C : Matrix (Fin 4) (Fin 4) ℝ → Bool) :
    ¬(C (distMatrix4 standardTet) = false ∧
      C (distMatrix4 (reflectZConfig standardTet)) = true) := by
  rw [distMatrix4_standardTet_reflect]
  intro ⟨h1, h2⟩
  simp_all

theorem no_distance_only_classifier_separates_reflection
    (P : Config4) (C : Matrix (Fin 4) (Fin 4) ℝ → Bool) :
    C (distMatrix4 (reflectZConfig P)) = C (distMatrix4 P) := by
  rw [distMatrix4_reflectZConfig]

theorem distance_only_function_reflection_invariant
    {α : Type} (F : Matrix (Fin 4) (Fin 4) ℝ → α) (P : Config4) :
    F (distMatrix4 (reflectZConfig P)) = F (distMatrix4 P) := by
  rw [distMatrix4_reflectZConfig]

-- ============================================================
-- SECTION 6. Distance-only score/loss impossibility
-- ============================================================

theorem distance_only_score_reflection_invariant
    (S : Matrix (Fin 4) (Fin 4) ℝ → ℝ) (P : Config4) :
    S (distMatrix4 (reflectZConfig P)) = S (distMatrix4 P) := by
  rw [distMatrix4_reflectZConfig]

theorem distance_only_score_cannot_penalize_mirror_standardTet
    (S : Matrix (Fin 4) (Fin 4) ℝ → ℝ) :
    S (distMatrix4 (reflectZConfig standardTet)) = S (distMatrix4 standardTet) := by
  rw [distMatrix4_standardTet_reflect]

-- ============================================================
-- SECTION 7. Combined witness theorem
-- ============================================================

theorem same_distances_opposite_orientations_standardTet :
    distMatrix4 (reflectZConfig standardTet) = distMatrix4 standardTet ∧
    orient4 standardTet = 1 ∧
    orient4 (reflectZConfig standardTet) = -1 :=
  ⟨distMatrix4_standardTet_reflect, orient4_standardTet, orient4_reflectZ_standardTet⟩

-- ============================================================
-- SECTION 8. Generalization to n points
-- ============================================================

abbrev ConfigN (n : ℕ) := Fin n → Point3

def distMatrixN (n : ℕ) (P : ConfigN n) : Matrix (Fin n) (Fin n) ℝ :=
  fun i j => distSq (P i) (P j)

def reflectZConfigN (n : ℕ) (P : ConfigN n) : ConfigN n :=
  fun k => reflectZPoint (P k)

theorem distMatrixN_reflectZConfigN (n : ℕ) (P : ConfigN n) :
    distMatrixN n (reflectZConfigN n P) = distMatrixN n P := by
  ext i j
  exact distSq_reflectZPoint (P i) (P j)

theorem distance_only_function_reflection_invariant_N
    {n : ℕ} {α : Type} (F : Matrix (Fin n) (Fin n) ℝ → α) (P : ConfigN n) :
    F (distMatrixN n (reflectZConfigN n P)) = F (distMatrixN n P) := by
  rw [distMatrixN_reflectZConfigN]