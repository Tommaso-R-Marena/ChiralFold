"""Regression tests for the correctness fixes shipped in v3.6.0.

Each test names the defect it locks down and asserts the *corrected* behaviour,
so a future refactor that reintroduces the bug fails here rather than silently
changing published numbers.
"""

from __future__ import annotations

import math
import textwrap

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chiralfold import _pdbio
from chiralfold import auditor as auditor_mod
from chiralfold import interface_scorer as iface_mod
from chiralfold import rotamers as rot_mod
from chiralfold.af3_correct import detect_chirality_violations
from chiralfold.enumerate import _random_patterns
from chiralfold.model import d_peptide_smiles, l_peptide_smiles, mixed_peptide_smiles
from chiralfold.pdb_pipeline import mirror_pdb
from chiralfold.structure_io import mmcif_to_pdb


# ─────────────────────────────────────────────────────────────────────────────
# model.py — C-terminal glycine emitted a spurious extra residue
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("length", [1, 2, 3, 4, 5])
def test_c_terminal_glycine_residue_count(length):
    """Poly-glycine must have exactly 4 heavy atoms per residue plus one OH.

    Up to v3.5.1 a C-terminal glycine appended two fragments — a malformed
    carbamate ``NC(=O)O`` *and* the real ``NCC(=O)O`` — so every sequence ending
    in G was built one residue too long (``d_peptide_smiles("G")`` produced a
    two-unit molecule).
    """
    seq = "G" * length
    mol = Chem.MolFromSmiles(d_peptide_smiles(seq))
    assert mol is not None
    assert mol.GetNumAtoms() == 4 * length + 1


def test_single_glycine_is_glycine():
    mol = Chem.MolFromSmiles(d_peptide_smiles("G"))
    assert rdMolDescriptors.CalcMolFormula(mol) == "C2H5NO2"


def test_glycine_anywhere_matches_alanine_backbone_length():
    """Swapping A→G must remove exactly one carbon, wherever G sits."""
    for probe in ("AAG", "AGA", "GAA"):
        gly = Chem.MolFromSmiles(d_peptide_smiles(probe))
        ala = Chem.MolFromSmiles(d_peptide_smiles("AAA"))
        assert gly.GetNumAtoms() == ala.GetNumAtoms() - 1, probe


def test_all_glycine_peptides_are_achiral():
    """Glycine has no Cα stereocentre, so the L and D builds must coincide."""
    for seq in ("G", "GG", "GGG", "GGGG"):
        assert Chem.CanonSmiles(d_peptide_smiles(seq)) == Chem.CanonSmiles(
            l_peptide_smiles(seq)
        )


def test_glycine_position_does_not_affect_other_residues_chirality():
    """Only the non-glycine residues may differ between the L and D builds."""
    for seq in ("AG", "GAG", "AGA"):
        d_mol = Chem.MolFromSmiles(d_peptide_smiles(seq))
        l_mol = Chem.MolFromSmiles(l_peptide_smiles(seq))
        # Same constitution, opposite configuration at every non-Gly residue.
        assert Chem.MolToSmiles(d_mol) != Chem.MolToSmiles(l_mol), seq
        assert Chem.MolToSmiles(
            Chem.MolFromSmiles(Chem.MolToSmiles(d_mol, isomericSmiles=False))
        ) == Chem.MolToSmiles(
            Chem.MolFromSmiles(Chem.MolToSmiles(l_mol, isomericSmiles=False))
        ), seq
        n_stereo = sum(1 for aa in seq if aa != "G")
        assert len(Chem.FindMolChiralCenters(d_mol)) == n_stereo, seq


def test_mixed_pattern_glycine_ignores_chirality_code():
    assert mixed_peptide_smiles("AG", "DD") == mixed_peptide_smiles("AG", "DL")


# ─────────────────────────────────────────────────────────────────────────────
# rotamers.py — chi1 sign convention was inverted
# ─────────────────────────────────────────────────────────────────────────────

def test_compute_chi1_matches_iupac_dihedral():
    """chi1 must equal the IUPAC/BioPython dihedral, not its negation.

    Up to v3.5.1 ``compute_chi1`` returned ``-dihedral``, so L side chains in
    the common g- rotamer were read as g+ and counted as outliers.
    """
    rng = np.random.default_rng(11)
    for _ in range(500):
        p = rng.normal(size=(4, 3))
        chi = rot_mod.compute_chi1(*p)
        ref = auditor_mod._dihedral_deg(*p)
        assert math.isclose(chi, ref, abs_tol=1e-9)


def test_chi1_gauche_minus_is_favored_for_l_residues():
    """A textbook g- leucine (chi1 ≈ -60°) must score favored, not outlier."""
    assert rot_mod._classify_chi1(-60.0, rot_mod.ROTAMER_LIBRARY["LEU"]) == "favored"
    assert rot_mod._classify_chi1(+60.0, rot_mod.ROTAMER_LIBRARY["DLE"]) == "favored"


def test_d_residues_have_mirrored_rotamer_wells():
    """D-CCD codes must be present and mirror their L counterparts."""
    for l3, d3 in (("LEU", "DLE"), ("PHE", "DPN"), ("LYS", "DLY"), ("MET", "MED")):
        assert d3 in rot_mod.ROTAMER_LIBRARY, d3
        l_wells = sorted(rot_mod.ROTAMER_LIBRARY[l3])
        d_wells = sorted(rot_mod.ROTAMER_LIBRARY[d3])
        mirrored = sorted(
            180.0 if abs(abs(w) - 180.0) < 1e-9 else -w for w in l_wells
        )
        assert d_wells == mirrored
        assert d3 in rot_mod._CG_ATOM_NAME


def test_threonine_gauche_plus_is_favored():
    """THR g+ is its most populated well (Lovell 2000) and must not be an outlier."""
    assert rot_mod._classify_chi1(62.0, rot_mod.ROTAMER_LIBRARY["THR"]) == "favored"


def test_validate_rotamers_scores_d_peptides(tmp_path):
    """A pure D structure must produce a non-empty rotamer report."""
    pdb = tmp_path / "d_leu.pdb"
    pdb.write_text(
        textwrap.dedent(
            """\
            HETATM    1  N   DLE A   1       1.320   0.000   0.000  1.00  0.00           N
            HETATM    2  CA  DLE A   1       0.000   0.000   0.000  1.00  0.00           C
            HETATM    3  CB  DLE A   1      -0.540   1.430   0.000  1.00  0.00           C
            HETATM    4  CG  DLE A   1      -0.100   2.200   1.260  1.00  0.00           C
            END
            """
        )
    )
    report = rot_mod.validate_rotamers(str(pdb))
    assert report["n_residues_checked"] == 1


# ─────────────────────────────────────────────────────────────────────────────
# _pdbio.py — altloc, insertion-code and multi-model handling
# ─────────────────────────────────────────────────────────────────────────────

def test_residue_modelled_only_as_altloc_b_is_kept(tmp_path):
    """A residue with no blank/'A' alternate must not vanish.

    The old accept-list ``{' ', 'A', '1'}`` dropped every atom of such a
    residue, silently deleting it from the audit.
    """
    pdb = tmp_path / "altloc_b_only.pdb"
    pdb.write_text(
        "ATOM      1  N  BALA A   1       0.000   0.000   0.000  0.60  0.00           N\n"
        "ATOM      2  CA CALA A   1       1.458   0.000   0.000  0.40  0.00           C\n"
        "END\n"
    )
    records = _pdbio.read_atom_records(str(pdb))
    assert len(records) == 1
    # Highest summed occupancy wins: B (0.60) beats C (0.40).
    assert records[0].altloc == "B"


def test_altloc_choice_is_residue_consistent(tmp_path):
    """All kept atoms of a residue must come from the same alternate."""
    pdb = tmp_path / "altloc_mixed.pdb"
    pdb.write_text(
        "ATOM      1  N  AALA A   1       0.000   0.000   0.000  0.40  0.00           N\n"
        "ATOM      2  N  BALA A   1       0.100   0.000   0.000  0.60  0.00           N\n"
        "ATOM      3  CA AALA A   1       1.458   0.000   0.000  0.40  0.00           C\n"
        "ATOM      4  CA BALA A   1       1.558   0.000   0.000  0.60  0.00           C\n"
        "END\n"
    )
    records = _pdbio.read_atom_records(str(pdb))
    assert {r.altloc for r in records} == {"B"}
    assert len(records) == 2


def test_insertion_codes_are_distinct_residues(tmp_path):
    """Residues 47 and 47A must not overwrite each other."""
    pdb = tmp_path / "icode.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A  47       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  CA  ALA A  47A      3.800   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    records = _pdbio.read_atom_records(str(pdb))
    assert len(records) == 2
    assert {r.icode for r in records} == {" ", "A"}
    # interface_scorer used to key de-duplication on (chain, resseq, name)
    assert len(iface_mod._parse_atoms(str(pdb))) == 2


def test_altloc_policy_only_affects_alternate_bearing_residues():
    """Structures with no lettered alternates must audit identically.

    This bounds the effect of the v3.6.0 per-residue altloc policy on any
    previously published aggregate: only residues that actually carry alternate
    conformations can change. ``benchmarks/altloc_policy_sensitivity.py``
    reports the measured deltas across the repository corpus.
    """
    from chiralfold.auditor import audit_pdb

    path = "chiralfold/data/examples/toy_ubiquitin_fragment.pdb"
    records = _pdbio.read_atom_records(path, resolve_altlocs=False)
    assert not [r for r in records if r.altloc not in _pdbio.BLANK_ALTLOCS]

    original = auditor_mod.read_atom_records

    def _legacy(p, **kwargs):
        kwargs["altloc_policy"] = "legacy"
        return original(p, **kwargs)

    auditor_mod.read_atom_records = _legacy
    try:
        legacy_report = audit_pdb(path)
    finally:
        auditor_mod.read_atom_records = original

    assert legacy_report == audit_pdb(path)


def test_unknown_altloc_policy_is_rejected():
    with pytest.raises(ValueError, match="unknown altloc policy"):
        _pdbio.select_altlocs([], policy="best-guess")


def test_only_first_model_is_read_by_default(tmp_path):
    pdb = tmp_path / "two_models.pdb"
    pdb.write_text(
        textwrap.dedent(
            """\
            MODEL        1
            ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
            ENDMDL
            MODEL        2
            ATOM      1  CA  ALA A   1      10.000   0.000   0.000  1.00  0.00           C
            ENDMDL
            END
            """
        )
    )
    first = _pdbio.read_atom_records(str(pdb))
    assert len(first) == 1
    assert first[0].x == pytest.approx(0.0)

    both = _pdbio.read_atom_records(str(pdb), first_model_only=False)
    assert len(both) == 2
    assert sorted(r.model for r in both) == [1, 2]


# ─────────────────────────────────────────────────────────────────────────────
# af3_correct.py — pseudo-Cβ, insertion codes, multi-model
# ─────────────────────────────────────────────────────────────────────────────

def _res_lines(resname, chain, resseq, icode=" ", cb_z=1.2, serial=1):
    return (
        f"ATOM  {serial:5d}  N   {resname} {chain}{resseq:4d}{icode}"
        "      1.201   0.847   0.000  1.00  0.00           N\n"
        f"ATOM  {serial + 1:5d}  CA  {resname} {chain}{resseq:4d}{icode}"
        "      0.000   0.000   0.000  1.00  0.00           C\n"
        f"ATOM  {serial + 2:5d}  C   {resname} {chain}{resseq:4d}{icode}"
        "     -1.250   0.881   0.000  1.00  0.00           C\n"
        f"ATOM  {serial + 3:5d}  CB  {resname} {chain}{resseq:4d}{icode}"
        f"      0.000  -0.500{cb_z:8.3f}  1.00  0.00           C\n"
    )


def test_d_residue_without_cb_is_not_a_false_violation(tmp_path):
    """A D-residue whose side chain is unmodelled must be unassignable.

    The old idealised-Cβ fallback made the signed volume positive by
    construction, so every such D-residue was reported as an L-coordinate
    violation.
    """
    pdb = tmp_path / "d_no_cb.pdb"
    pdb.write_text(
        "HETATM    1  N   DAL A   1       1.201   0.847   0.000  1.00  0.00           N\n"
        "HETATM    2  CA  DAL A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    3  C   DAL A   1      -1.250   0.881   0.000  1.00  0.00           C\n"
        "END\n"
    )
    report = detect_chirality_violations(str(pdb))
    assert report["n_violations"] == 0
    assert report["n_unassignable"] == 1


def test_insertion_code_residues_are_classified_independently(tmp_path):
    """A violation at residue 1 must not be attributed to residue 1A too."""
    pdb = tmp_path / "icode_chirality.pdb"
    pdb.write_text(
        _res_lines("ALA", "A", 1, " ", cb_z=-1.2, serial=1)      # inverted → violation
        + _res_lines("ALA", "A", 1, "A", cb_z=+1.2, serial=11)   # correct L
        + "END\n"
    )
    report = detect_chirality_violations(str(pdb))
    assert report["n_checked"] == 2
    assert report["n_violations"] == 1
    assert report["violations"][0]["icode"] == " "


def test_all_models_are_detected_and_corrected(tmp_path):
    """Violations in models 2..n must be found and fixed, not conflated."""
    body = _res_lines("ALA", "A", 1, " ", cb_z=-1.2, serial=1)
    pdb = tmp_path / "multi_model.pdb"
    pdb.write_text(
        "MODEL        1\n" + body + "ENDMDL\n"
        "MODEL        2\n" + body + "ENDMDL\nEND\n"
    )
    report = detect_chirality_violations(str(pdb))
    assert report["n_models"] == 2
    assert report["n_violations"] == 2
    assert sorted(v["model"] for v in report["violations"]) == [1, 2]

    from chiralfold.af3_correct import correct_af3_output

    out = tmp_path / "fixed.pdb"
    result = correct_af3_output(str(pdb), str(out))
    assert result["correction"]["n_corrected"] == 2
    assert result["after"]["n_violations"] == 0


# ─────────────────────────────────────────────────────────────────────────────
# pdb_pipeline.py — mirroring an NMR ensemble flattened the models
# ─────────────────────────────────────────────────────────────────────────────

def test_mirror_preserves_model_records(tmp_path):
    """An n-model input must produce an n-model output, not n stacked copies."""
    body = (
        "ATOM      1  N   ALA A   1       1.201   0.847   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   ALA A   1      -1.250   0.881   0.000  1.00  0.00           C\n"
    )
    src = tmp_path / "ensemble.pdb"
    src.write_text(
        "MODEL        1\n" + body + "ENDMDL\n"
        "MODEL        2\n" + body + "ENDMDL\nEND\n"
    )
    out = tmp_path / "ensemble_D.pdb"
    result = mirror_pdb(str(src), str(out))

    assert result["n_models"] == 2
    assert result["n_residues"] == 1          # per model, not 2
    text = out.read_text()
    assert text.count("MODEL ") == 2
    assert text.count("ENDMDL") == 2

    # Downstream readers must therefore see one residue, not two copies.
    records = _pdbio.read_atom_records(str(out))
    assert len({r.residue_key for r in records}) == 1


def test_mirror_keeps_cryst_and_title_records(tmp_path):
    src = tmp_path / "with_header.pdb"
    src.write_text(
        "TITLE     TEST STRUCTURE\n"
        "CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1           1\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    out = tmp_path / "mirrored.pdb"
    mirror_pdb(str(src), str(out))
    text = out.read_text()
    assert "TITLE     TEST STRUCTURE" in text
    assert "CRYST1" in text


# ─────────────────────────────────────────────────────────────────────────────
# interface_scorer.py — salt bridges were measured Cα–Cα at 4 Å
# ─────────────────────────────────────────────────────────────────────────────

def test_salt_bridge_uses_charged_side_chain_groups(tmp_path):
    """LYS NZ ··· ASP OD1 at 3.2 Å is a salt bridge; Cα–Cα is 9 Å apart.

    With the old Cα–Cα 4 Å rule this real ion pair scored zero, because Cα atoms
    across an interface are essentially never within 4 Å.
    """
    receptor = tmp_path / "rec.pdb"
    ligand = tmp_path / "lig.pdb"
    receptor.write_text(
        "ATOM      1  CA  LYS A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  NZ  LYS A   1       4.000   0.000   0.000  1.00  0.00           N\n"
        "END\n"
    )
    ligand.write_text(
        "ATOM      1  CA  ASP B   1       9.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  OD1 ASP B   1       7.200   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    rec = iface_mod._parse_atoms(str(receptor))
    lig = iface_mod._parse_atoms(str(ligand))
    assert iface_mod._compute_salt_bridges(rec, lig) == 1


def test_salt_bridge_counted_once_per_residue_pair(tmp_path):
    """Arg has three cationic N atoms; one ion pair must still count once."""
    receptor = tmp_path / "arg.pdb"
    ligand = tmp_path / "glu.pdb"
    receptor.write_text(
        "ATOM      1  NE  ARG A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  NH1 ARG A   1       0.500   0.500   0.000  1.00  0.00           N\n"
        "ATOM      3  NH2 ARG A   1       0.500  -0.500   0.000  1.00  0.00           N\n"
        "END\n"
    )
    ligand.write_text(
        "ATOM      1  OE1 GLU B   1       2.800   0.000   0.000  1.00  0.00           O\n"
        "ATOM      2  OE2 GLU B   1       3.000   0.800   0.000  1.00  0.00           O\n"
        "END\n"
    )
    rec = iface_mod._parse_atoms(str(receptor))
    lig = iface_mod._parse_atoms(str(ligand))
    assert iface_mod._compute_salt_bridges(rec, lig) == 1


def test_far_charged_groups_are_not_a_salt_bridge(tmp_path):
    receptor = tmp_path / "rec_far.pdb"
    ligand = tmp_path / "lig_far.pdb"
    receptor.write_text(
        "ATOM      1  NZ  LYS A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n"
    )
    ligand.write_text(
        "ATOM      1  OD1 ASP B   1      12.000   0.000   0.000  1.00  0.00           O\nEND\n"
    )
    assert iface_mod._compute_salt_bridges(
        iface_mod._parse_atoms(str(receptor)),
        iface_mod._parse_atoms(str(ligand)),
    ) == 0


def test_interface_pairs_match_brute_force(tmp_path):
    """The KD-tree contact search must reproduce the dense-matrix result."""
    rng = np.random.default_rng(5)
    rec_xyz = rng.uniform(0, 12, size=(60, 3))
    lig_xyz = rng.uniform(6, 18, size=(50, 3))

    def _write(path, xyz, chain):
        lines = [
            f"ATOM  {i + 1:5d}  CA  ALA {chain}{i + 1:4d}    "
            f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00           C\n"
            for i, p in enumerate(xyz)
        ]
        path.write_text("".join(lines) + "END\n")

    rec_path = tmp_path / "rec_rand.pdb"
    lig_path = tmp_path / "lig_rand.pdb"
    _write(rec_path, rec_xyz, "A")
    _write(lig_path, lig_xyz, "B")

    rec = iface_mod._parse_atoms(str(rec_path))
    lig = iface_mod._parse_atoms(str(lig_path))
    pairs = iface_mod._find_interface_pairs(rec, lig, cutoff=5.0)

    dense = np.linalg.norm(rec_xyz[:, None, :] - lig_xyz[None, :, :], axis=-1)
    expected = int((dense <= 5.0).sum())
    assert len(pairs) == expected
    for _r, _l, d in pairs:
        assert d <= 5.0 + 1e-9


# ─────────────────────────────────────────────────────────────────────────────
# enumerate.py — rejection sampling could return too few patterns
# ─────────────────────────────────────────────────────────────────────────────

def test_random_patterns_returns_exact_count():
    import random

    rng = random.Random(0)
    patterns = _random_patterns(n_stereo=12, n_samples=200, rng=rng)
    assert len(patterns) == 200
    assert len(set(patterns)) == 200
    assert all(len(p) == 12 for p in patterns)
    assert all(set(p) <= {"L", "D"} for p in patterns)


def test_random_patterns_exhaustive_when_space_is_small():
    import random

    patterns = _random_patterns(n_stereo=3, n_samples=100, rng=random.Random(0))
    assert len(patterns) == 8
    assert len(set(patterns)) == 8


# ─────────────────────────────────────────────────────────────────────────────
# structure_io.py — quoted mmCIF atom names kept their quotes
# ─────────────────────────────────────────────────────────────────────────────

def test_mmcif_quoted_atom_names_are_unquoted():
    cif = textwrap.dedent(
        """\
        data_test
        loop_
        _atom_site.group_PDB
        _atom_site.id
        _atom_site.type_symbol
        _atom_site.label_atom_id
        _atom_site.label_comp_id
        _atom_site.label_asym_id
        _atom_site.label_seq_id
        _atom_site.Cartn_x
        _atom_site.Cartn_y
        _atom_site.Cartn_z
        _atom_site.occupancy
        _atom_site.B_iso_or_equiv
        _atom_site.pdbx_PDB_model_num
        ATOM 1 O "O5'" DA A 1 1.000 2.000 3.000 1.00 20.00 1
        ATOM 2 C "C1'" DA A 1 2.000 2.000 3.000 1.00 20.00 1
        """
    )
    pdb = mmcif_to_pdb(cif)
    assert pdb is not None
    names = [line[12:16].strip() for line in pdb.splitlines() if line.startswith("ATOM")]
    assert names == ["O5'", "C1'"]
    assert not any('"' in n for n in names)


# ─────────────────────────────────────────────────────────────────────────────
# auditor.py — vectorised kernels must match the scalar reference
# ─────────────────────────────────────────────────────────────────────────────

def test_batched_geometry_kernels_match_scalar():
    rng = np.random.default_rng(17)
    pts = rng.normal(size=(200, 4, 3))

    dih_v = auditor_mod._dihedrals_v(
        pts[:, 0], pts[:, 1], pts[:, 2], pts[:, 3]
    )
    ang_v = auditor_mod._bond_angles_v(pts[:, 0], pts[:, 1], pts[:, 2])
    len_v = auditor_mod._bond_lengths_v(pts[:, 0], pts[:, 1])
    vol_v = auditor_mod._signed_volumes_v(
        pts[:, 0], pts[:, 1], pts[:, 2], pts[:, 3]
    )

    for i, p in enumerate(pts):
        assert dih_v[i] == pytest.approx(auditor_mod._dihedral_deg(*p), abs=1e-9)
        assert ang_v[i] == pytest.approx(
            auditor_mod._bond_angle_deg(p[0], p[1], p[2]), abs=1e-9
        )
        assert len_v[i] == pytest.approx(
            auditor_mod._bond_length(p[0], p[1]), abs=1e-12
        )
        assert vol_v[i] == pytest.approx(auditor_mod._signed_volume(*p), abs=1e-12)


def test_batched_hydrogen_placement_matches_scalar():
    rng = np.random.default_rng(23)
    ca = rng.normal(size=(40, 3))
    n = ca + rng.normal(size=(40, 3))
    c = ca + rng.normal(size=(40, 3))
    cb = ca + rng.normal(size=(40, 3))

    batch = auditor_mod._place_tetrahedral_h_batch(ca, (n, c, cb), 1.09)
    for i in range(ca.shape[0]):
        ref = auditor_mod._place_tetrahedral_h(ca[i], [n[i], c[i], cb[i]], 1.09)
        assert np.allclose(batch[i], ref, atol=1e-12)

    ha2, ha3 = auditor_mod._place_gly_ha2_ha3_batch(n, ca, c)
    for i in range(ca.shape[0]):
        r2, r3 = auditor_mod._place_gly_ha2_ha3(n[i], ca[i], c[i])
        assert np.allclose(ha2[i], r2, atol=1e-12)
        assert np.allclose(ha3[i], r3, atol=1e-12)


def test_classify_regions_matches_scalar_scorer():
    from chiralfold.ramachandran import classify_regions, score_ramachandran

    rng = np.random.default_rng(31)
    phi = rng.uniform(-180, 180, size=400)
    psi = rng.uniform(-180, 180, size=400)
    rtypes = list(
        rng.choice(
            ["general", "glycine", "proline", "D-general", "D-proline"], size=400
        )
    )
    batch = classify_regions(phi, psi, rtypes)
    for i in range(len(batch)):
        assert batch[i] == score_ramachandran(phi[i], psi[i], rtypes[i])


def test_exclusion_keys_and_set_agree():
    atoms = auditor_mod._parse_pdb("chiralfold/data/examples/toy_ubiquitin_fragment.pdb")
    heavy = [a for a in atoms if not auditor_mod._is_hydrogen_atom(a)]
    with_h = auditor_mod._add_backbone_hydrogens(heavy)
    n = len(with_h)

    keys = auditor_mod._excluded_pair_keys(with_h)
    as_set = auditor_mod._clash_excluded_index_pairs(with_h)

    assert len(keys) == len(as_set)
    assert np.all(np.diff(keys) > 0)  # sorted and unique
    assert {(int(k // n), int(k % n)) for k in keys} == as_set


def test_clash_query_radius_cannot_miss_a_clash():
    """The tightened KD-tree radius must still cover every possible overlap."""
    max_r = max(auditor_mod.VDW_RADII.values())
    assert max_r >= auditor_mod.VDW_DEFAULT
    widest_clash = 2 * max_r - auditor_mod.CLASH_OVERLAP_CUTOFF
    for r1 in auditor_mod.VDW_RADII.values():
        for r2 in auditor_mod.VDW_RADII.values():
            assert r1 + r2 - auditor_mod.CLASH_OVERLAP_CUTOFF <= widest_clash + 1e-12
