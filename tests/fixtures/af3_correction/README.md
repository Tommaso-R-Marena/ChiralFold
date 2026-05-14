# AF3 Correction Test Fixtures

These fixtures are **synthetic** PDB snippets generated programmatically by
`tests/test_af3_correction.py`. They are not derived from real AF3 outputs
because AlphaFold 3 weights are not publicly downloadable without an AlphaFold
Server account, and the Childs et al. (2025) D-peptide predictions are not
redistributable.

Two minimal residues are used as test inputs, each with deliberate chirality
inversions encoded by flipping the sign of the Cβ z-coordinate:

- `synthetic_l_ala_inverted.pdb` — ATOM record (L-Ala by label) with the
  Cβ placed on the wrong side of the N–Cα–C plane, producing a positive
  signed volume that the auditor flags as a chirality violation.
- `synthetic_d_ala_inverted.pdb` — HETATM record (DAL by label) with the Cβ
  placed on the L side of the N–Cα–C plane (negative signed volume), again a
  violation against the expected D-chirality.
- `synthetic_correct_l_ala.pdb` — a sanity-control L-Ala with correct chirality
  (negative signed volume), used to assert that `correct_af3_output` is a
  no-op on already-correct inputs.

If real AF3 predictions become available without account credentials in a
future revision, they can be dropped into this directory and the test can be
extended via the `@pytest.mark.parametrize` "real_af3_predictions" list.
