/-
# Concrete tetrahedron and mirror examples

These examples instantiate the general results with `standardTet` and normalize
the resulting numerical orientations.
-/
import Mathlib
import RequestProject.ModelAbstraction

/-- The concrete tetrahedron with vertices `0, e₁, e₂, e₃`. -/
def standardTetExample : Config4 := standardTet

/-- Its mirror image across the xy-plane. -/
def mirroredStandardTetExample : Config4 := reflectZConfig standardTetExample

/-- The concrete example and its mirror have equal pairwise distance matrices. -/
theorem standardTetExample_mirror_distances :
    distMatrix4 mirroredStandardTetExample = distMatrix4 standardTetExample := by
  norm_num [mirroredStandardTetExample, distMatrix4_reflectZConfig]

/-- Normalization gives positive orientation for the example. -/
theorem standardTetExample_orientation :
    orient4 standardTetExample = 1 := by
  norm_num [standardTetExample, orient4_standardTet]

/-- Normalization gives negative orientation for the mirror. -/
theorem mirroredStandardTetExample_orientation :
    orient4 mirroredStandardTetExample = -1 := by
  norm_num [mirroredStandardTetExample, standardTetExample,
    orient4_reflectZConfig, orient4_standardTet]

/-- Thus the concrete pair has matching distance matrices and differing signed
orientations. -/
theorem standardTetExample_same_distances_different_orientations :
    distMatrix4 mirroredStandardTetExample = distMatrix4 standardTetExample ∧
      orient4 mirroredStandardTetExample ≠ orient4 standardTetExample := by
  constructor
  · exact standardTetExample_mirror_distances
  · rw [mirroredStandardTetExample_orientation, standardTetExample_orientation]
    norm_num
