"""
Tests for chiralfold.af3_correct — AlphaFold 3 chirality correction pipeline.

We use synthetic PDB fixtures with deliberate chirality inversions because AF3
weights are not publicly redistributable. The fixture builders are isolated in
``tests/fixtures/af3_correction/`` and documented in the directory README.

For each fixture we check three properties of the correction pipeline:
  1. violations are detected before correction;
  2. the corrected coordinates produce a geometrically valid Cα (CA-N, CA-C,
     and CA-Cβ bond lengths within tolerance of the input);
  3. no new violations are introduced (in particular at neighbouring residues).
"""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import numpy as np
import pytest

from chiralfold.af3_correct import (
    detect_chirality_violations,
    correct_af3_output,
    correct_chirality,
    _estimate_cb,
    _format_coord,
    _parse_pdb_full,
    _reflect_across_plane,
    _rewrite_atom_coords,
    _signed_volume,
)

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "af3_correction"
FIXTURE_DIR.mkdir(parents=True, exist_ok=True)


# --- Synthetic PDB builders (deterministic; deliberate chirality flips) ----

_L_ALA_INVERTED = textwrap.dedent("""\
    ATOM      1  N   ALA A   1       1.201   0.847   0.000  1.00  0.00           N
    ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
    ATOM      3  C   ALA A   1      -1.250   0.881   0.000  1.00  0.00           C
    ATOM      4  O   ALA A   1      -1.200   2.095   0.000  1.00  0.00           O
    ATOM      5  CB  ALA A   1       0.000  -0.500  -1.200  1.00  0.00           C
    ATOM      6  HA  ALA A   1       0.100   0.100   1.000  1.00  0.00           H
    END
""")

_L_ALA_CORRECT = textwrap.dedent("""\
    ATOM      1  N   ALA A   1       1.201   0.847   0.000  1.00  0.00           N
    ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
    ATOM      3  C   ALA A   1      -1.250   0.881   0.000  1.00  0.00           C
    ATOM      4  O   ALA A   1      -1.200   2.095   0.000  1.00  0.00           O
    ATOM      5  CB  ALA A   1       0.000  -0.500   1.200  1.00  0.00           C
    END
""")

# DAL HETATM with Cβ on the L (negative-z) side: signed volume should come out
# positive -> "L found" against "D expected" -> a chirality violation.
_D_ALA_INVERTED = textwrap.dedent("""\
    HETATM    1  N   DAL B   1       1.201   0.847   0.000  1.00  0.00           N
    HETATM    2  CA  DAL B   1       0.000   0.000   0.000  1.00  0.00           C
    HETATM    3  C   DAL B   1      -1.250   0.881   0.000  1.00  0.00           C
    HETATM    4  O   DAL B   1      -1.200   2.095   0.000  1.00  0.00           O
    HETATM    5  CB  DAL B   1       0.000  -0.500   1.200  1.00  0.00           C
    END
""")


@pytest.fixture(scope="module")
def fixtures(tmp_path_factory) -> dict:
    """Materialize synthetic PDBs into both the repo fixture dir (for ground-
    truth artifacts) and a tmp scratch dir (for output files)."""
    scratch = tmp_path_factory.mktemp("af3_corr_outputs")

    layout = {
        "l_ala_inverted":  _L_ALA_INVERTED,
        "l_ala_correct":   _L_ALA_CORRECT,
        "d_ala_inverted":  _D_ALA_INVERTED,
    }
    paths = {}
    for stem, content in layout.items():
        fixture_path = FIXTURE_DIR / f"synthetic_{stem}.pdb"
        fixture_path.write_text(content)
        paths[stem] = fixture_path

    return {"fixtures": paths, "scratch": scratch}


# ----- Detection ----------------------------------------------------------

def test_detects_l_ala_inverted(fixtures):
    report = detect_chirality_violations(str(fixtures["fixtures"]["l_ala_inverted"]))
    assert report["n_checked"] == 1
    assert report["n_violations"] == 1
    v = report["violations"][0]
    assert v["resname"] == "ALA"
    assert v["expected_chirality"] == "L"
    assert v["found_chirality"] == "D"


def test_detects_d_ala_inverted(fixtures):
    report = detect_chirality_violations(str(fixtures["fixtures"]["d_ala_inverted"]))
    assert report["n_checked"] == 1
    assert report["n_violations"] == 1
    v = report["violations"][0]
    assert v["resname"] == "DAL"
    assert v["expected_chirality"] == "D"
    assert v["found_chirality"] == "L"


def test_clean_l_ala_has_no_violations(fixtures):
    report = detect_chirality_violations(str(fixtures["fixtures"]["l_ala_correct"]))
    assert report["n_violations"] == 0


# ----- Correction --------------------------------------------------------

@pytest.mark.parametrize("stem", ["l_ala_inverted", "d_ala_inverted"])
def test_correction_resolves_violations(fixtures, stem):
    inp = fixtures["fixtures"][stem]
    out = fixtures["scratch"] / f"{stem}_corrected.pdb"
    result = correct_af3_output(str(inp), str(out))

    assert result["before"]["n_violations"] >= 1
    assert result["after"]["n_violations"] == 0
    assert result["success"] is True
    assert result["correction"]["n_corrected"] >= 1
    assert out.exists()


def test_corrected_geometry_preserves_bond_lengths(fixtures):
    inp = fixtures["fixtures"]["l_ala_inverted"]
    out = fixtures["scratch"] / "l_ala_inverted_geom.pdb"
    correct_af3_output(str(inp), str(out))

    coords_before = _read_atom_xyz(inp)
    coords_after = _read_atom_xyz(out)

    # N, C, O, CB positions are not reflected — only CA (and HA if present) is.
    for name in ("N", "C", "O", "CB"):
        np.testing.assert_allclose(
            coords_after[name], coords_before[name], atol=1e-6,
            err_msg=f"non-CA atom {name} should be unchanged after correction"
        )

    # CA-N and CA-C distances should be preserved (reflection across the
    # plane defined by N, C, CB).
    for partner in ("N", "C", "CB"):
        d_before = np.linalg.norm(coords_before["CA"] - coords_before[partner])
        d_after = np.linalg.norm(coords_after["CA"] - coords_after[partner])
        assert abs(d_before - d_after) < 1e-3, (
            f"CA-{partner} bond length changed: {d_before:.4f} -> {d_after:.4f}"
        )

    # Signed volume should flip sign cleanly.
    sv_before = _signed_volume(
        coords_before["CA"], coords_before["N"],
        coords_before["C"], coords_before["CB"],
    )
    sv_after = _signed_volume(
        coords_after["CA"], coords_after["N"],
        coords_after["C"], coords_after["CB"],
    )
    assert sv_before * sv_after < 0, (
        f"Signed volume did not invert: {sv_before:.4f} -> {sv_after:.4f}"
    )


def test_correction_is_noop_on_clean_input(fixtures):
    inp = fixtures["fixtures"]["l_ala_correct"]
    out = fixtures["scratch"] / "clean_passthrough.pdb"
    result = correct_af3_output(str(inp), str(out))
    assert result["before"]["n_violations"] == 0
    assert result["after"]["n_violations"] == 0
    assert result["correction"]["n_corrected"] == 0


def test_parse_pdb_full_skips_water_altloc_duplicates_and_short_lines(tmp_path):
    pdb = tmp_path / "parse_edges.pdb"
    pdb.write_text(
        textwrap.dedent(
            """\
            ATOM      1  N   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
            ATOM      2  CA AALA A   2       0.000   0.000   0.000  1.00  0.00           C
            ATOM      3  CA BALA A   2       1.000   0.000   0.000  1.00  0.00           C
            ATOM      4  CA AALA A   2       2.000   0.000   0.000  1.00  0.00           C
            ATOM  broken
            END
            """
        )
    )

    lines, atoms = _parse_pdb_full(str(pdb))

    assert len(lines) == 6
    assert len(atoms) == 1
    assert atoms[0].altloc == "A"
    assert atoms[0].line_idx == 1


def test_detect_skips_gly_pro_missing_backbone_and_uses_estimated_cb(tmp_path):
    pdb = tmp_path / "detect_edges.pdb"
    pdb.write_text(
        textwrap.dedent(
            """\
            ATOM      1  N   GLY A   1       1.201   0.847   0.000  1.00  0.00           N
            ATOM      2  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C
            ATOM      3  C   GLY A   1      -1.250   0.881   0.000  1.00  0.00           C
            ATOM      4  N   PRO A   2       1.201   0.847   0.000  1.00  0.00           N
            ATOM      5  CA  PRO A   2       0.000   0.000   0.000  1.00  0.00           C
            ATOM      6  C   PRO A   2      -1.250   0.881   0.000  1.00  0.00           C
            ATOM      7  CA  ALA A   3       0.000   0.000   0.000  1.00  0.00           C
            HETATM    8  N   UNK B   1       1.201   0.847   0.000  1.00  0.00           N
            HETATM    9  CA  UNK B   1       0.000   0.000   0.000  1.00  0.00           C
            HETATM   10  C   UNK B   1      -1.250   0.881   0.000  1.00  0.00           C
            END
            """
        )
    )

    report = detect_chirality_violations(str(pdb))

    assert report["n_residues"] == 4
    assert report["n_checked"] == 1


def test_detect_missing_file_and_correct_chirality_missing_file_raise(tmp_path):
    missing = tmp_path / "missing.pdb"
    with pytest.raises(FileNotFoundError):
        detect_chirality_violations(str(missing))
    with pytest.raises(FileNotFoundError):
        correct_chirality(str(missing), str(tmp_path / "out.pdb"))
    with pytest.raises(FileNotFoundError):
        correct_af3_output(str(missing), str(tmp_path / "out.pdb"))


def test_correct_chirality_noop_copies_clean_input(fixtures):
    inp = fixtures["fixtures"]["l_ala_correct"]
    out = fixtures["scratch"] / "direct_noop.pdb"

    result = correct_chirality(str(inp), str(out))

    assert result == {
        "n_corrected": 0,
        "violations_before": 0,
        "violations_after": 0,
        "output_path": str(out),
    }
    assert out.read_text() == Path(inp).read_text()


def test_correct_af3_output_derives_default_output_path(tmp_path, fixtures):
    inp = tmp_path / "clean_copy.pdb"
    shutil.copy2(fixtures["fixtures"]["l_ala_correct"], inp)

    result = correct_af3_output(str(inp))

    assert result["output_path"].endswith("_corrected.pdb")
    assert result["success"] is True
    assert Path(result["output_path"]).exists()


def test_af3_geometry_helpers_handle_edge_cases():
    line = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n"

    assert _format_coord(1.25) == "   1.250"
    assert _rewrite_atom_coords("short", 9.0, 8.0, 7.0) == "short"
    rewritten = _rewrite_atom_coords(line, -1.0, -2.0, -3.0)
    assert float(rewritten[30:38]) == -1.0
    assert float(rewritten[38:46]) == -2.0
    assert float(rewritten[46:54]) == -3.0

    point = np.array([1.0, 2.0, 3.0])
    reflected = _reflect_across_plane(
        point,
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([2.0, 0.0, 0.0]),
    )
    np.testing.assert_allclose(reflected, point)

    cb = _estimate_cb(
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([-1.0, 0.0, 0.0]),
    )
    assert cb.shape == (3,)


# ----- Helpers -----------------------------------------------------------

def _read_atom_xyz(path) -> dict:
    """Parse a tiny PDB and return a {atom_name: np.array([x, y, z])} map."""
    coords = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords[name] = np.array([x, y, z])
    return coords
