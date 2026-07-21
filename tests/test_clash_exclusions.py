"""Regression tests for MolProbity-like clash exclusions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chiralfold import auditor as auditor_mod
from chiralfold.auditor import audit_pdb


TOY = (
    Path(__file__).resolve().parents[1]
    / "chiralfold"
    / "data"
    / "examples"
    / "toy_ubiquitin_fragment.pdb"
)


def test_leu_ca_cg_is_excluded_as_1_3():
    """CA–CG (~2.6 Å) is covalent 1-3 via CB — must never score as a clash."""
    atoms = auditor_mod._parse_pdb(str(TOY))
    leus = {}
    for a in atoms:
        if a.resname == "LEU" and a.name.strip() in ("CA", "CG", "CB"):
            leus.setdefault(a.resseq, {})[a.name.strip()] = a
    assert leus, "toy fragment should contain LEU"
    for resseq, m in leus.items():
        if "CA" in m and "CG" in m:
            assert auditor_mod._are_bonded_or_angled(m["CA"], m["CG"]) is True


def test_amide_hydrogen_placement_not_on_prev_ca():
    atoms = auditor_mod._parse_pdb(str(TOY))
    with_h = auditor_mod._add_backbone_hydrogens(atoms)
    # No proline should receive a fake amide H
    for a in with_h:
        if a.element_upper == "H" and a.resname in auditor_mod.PROLINE_RESNAMES:
            pytest.fail(f"Proline {a.resseq} must not get an amide H")
    # Previous CA to amide H should be a normal 1-4 distance (~2.0–2.8 Å)
    by = {}
    for a in with_h:
        by.setdefault((a.chain, a.resseq), {})[a.name.strip()] = a
    keys = sorted(by)
    checked = 0
    for (c0, r0), (c1, r1) in zip(keys, keys[1:]):
        if c0 != c1 or r1 - r0 != 1:
            continue
        if "CA" not in by[(c0, r0)] or "H" not in by[(c1, r1)]:
            continue
        dist = float(
            np.linalg.norm(by[(c0, r0)]["CA"].xyz - by[(c1, r1)]["H"].xyz)
        )
        assert dist > 2.0, f"amide H too close to prev CA: {dist}"
        checked += 1
    assert checked >= 1


def test_toy_ubiquitin_clashscore_is_reasonable():
    """After topology exclusions, a clean fragment must not look clash-ridden."""
    report = audit_pdb(str(TOY))
    assert report["clashes"]["clash_score"] < 50.0
    assert report["chirality"]["n_wrong"] == 0


@pytest.mark.slow
def test_experimental_1ubq_clashscore_ballpark():
    from chiralfold.fetch import fetch_structure

    path = fetch_structure("1UBQ").path
    report = audit_pdb(path)
    # Pre-fix scores were >200 from false 1-3 hits; MolProbity is ~0–5.
    assert report["clashes"]["clash_score"] < 40.0
    assert report["planarity"]["pct_within_6deg"] > 90.0
    assert report["ramachandran"]["pct_favored"] > 90.0
