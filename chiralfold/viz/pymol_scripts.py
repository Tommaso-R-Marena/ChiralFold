"""
PyMOL script writers for ChiralFold structure visualization.

These helpers write plain-text ``.pml`` files that can be opened with PyMOL:

    pymol chirality_audit.pml
    pymol -cq chirality_audit.pml

PyMOL is intentionally not imported here. Install system PyMOL separately when
rendering sessions or PNGs; ChiralFold's package dependency set stays light.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from chiralfold.af3_correct import detect_chirality_violations

ReportLike = Union[Dict[str, object], Sequence[Dict[str, object]], str, os.PathLike]
ResidueKey = Tuple[str, int]


def _quote_path(path: Union[str, os.PathLike]) -> str:
    return str(Path(path)).replace("\\", "/")


def _object_name(path: Union[str, os.PathLike], suffix: str = "") -> str:
    stem = Path(path).stem
    safe = "".join(ch if ch.isalnum() else "_" for ch in stem)
    return f"{safe}{suffix}"


def _load_report(report_or_violations: ReportLike) -> Union[Dict[str, object], Sequence[Dict[str, object]]]:
    if isinstance(report_or_violations, (str, os.PathLike)):
        path = Path(report_or_violations)
        with path.open() as fh:
            return json.load(fh)
    return report_or_violations


def _extract_violations(report_or_violations: ReportLike) -> List[Dict[str, object]]:
    data = _load_report(report_or_violations)
    if isinstance(data, dict):
        if "violations" in data:
            return list(data.get("violations", []))
        if "chirality" in data and isinstance(data["chirality"], dict):
            return list(data["chirality"].get("violations", []))
        if "before" in data and isinstance(data["before"], dict):
            return list(data["before"].get("violations", []))
    return list(data)


def _violation_residues(report_or_violations: ReportLike) -> List[ResidueKey]:
    residues = []
    for violation in _extract_violations(report_or_violations):
        chain = str(violation.get("chain", "")).strip()
        resnum = violation.get("resnum", violation.get("resseq", None))
        if resnum is None and "residue" in violation:
            label = str(violation["residue"])
            digits = "".join(ch for ch in label if ch.isdigit())
            resnum = int(digits) if digits else None
            if ":" in label and not chain:
                chain = label.split(":", 1)[0].strip()
        if resnum is None:
            continue
        residues.append((chain, int(resnum)))
    return sorted(set(residues), key=lambda item: (item[0], item[1]))


def _selection_for_residues(object_name: str, residues: Iterable[ResidueKey]) -> str:
    parts = []
    for chain, resnum in residues:
        if chain:
            parts.append(f"({object_name} and chain {chain} and resi {resnum})")
        else:
            parts.append(f"({object_name} and resi {resnum})")
    return " or ".join(parts) if parts else "none"


def _write_lines(out_pml: Union[str, os.PathLike], lines: Sequence[str]) -> Path:
    path = Path(out_pml)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def write_chirality_audit_session(
    pdb_path: Union[str, os.PathLike],
    report_or_violations: ReportLike,
    out_pml: Union[str, os.PathLike],
) -> Path:
    """
    Write a PyMOL script coloring residues by chirality audit status.

    Correct/checkable residues are green, residues with chirality violations are
    red, and glycine is gray. ``report_or_violations`` may be a detector report,
    an auditor report, a list of violation dictionaries, or a JSON path.
    """
    obj = _object_name(pdb_path, "_audit")
    wrong_selection = _selection_for_residues(obj, _violation_residues(report_or_violations))
    lines = [
        "reinitialize",
        "bg_color white",
        f"load {_quote_path(pdb_path)}, {obj}",
        f"hide everything, {obj}",
        f"show cartoon, {obj}",
        f"show sticks, {obj} and not elem H",
        f"color forest, {obj}",
        f"color gray70, {obj} and resn GLY",
        f"select chirality_wrong, {wrong_selection}",
        "color red, chirality_wrong",
        "show sticks, chirality_wrong",
        "show spheres, chirality_wrong and name CA+CB",
        'label chirality_wrong and name CA, "%s%s" % (resn, resi)',
        "set sphere_scale, 0.35",
        "set cartoon_fancy_helices, 1",
        "set ray_opaque_background, off",
        f"zoom {obj}",
    ]
    return _write_lines(out_pml, lines)


def write_mirror_comparison_session(
    l_pdb: Union[str, os.PathLike],
    d_pdb: Union[str, os.PathLike],
    out_pml: Union[str, os.PathLike],
) -> Path:
    """Write a side-by-side PyMOL script comparing L and D structures."""
    l_obj = _object_name(l_pdb, "_L")
    d_obj = _object_name(d_pdb, "_D")
    lines = [
        "reinitialize",
        "bg_color white",
        f"load {_quote_path(l_pdb)}, {l_obj}",
        f"load {_quote_path(d_pdb)}, {d_obj}",
        f"hide everything, {l_obj} or {d_obj}",
        f"show cartoon, {l_obj} or {d_obj}",
        f"show sticks, ({l_obj} or {d_obj}) and not elem H",
        f"color marine, {l_obj}",
        f"color magenta, {d_obj}",
        f"translate [35, 0, 0], {d_obj}",
        f"label {l_obj} and name CA and resi 1, 'L structure'",
        f"label {d_obj} and name CA and resi 1, 'D/mirror structure'",
        "set cartoon_transparency, 0.08",
        "set ray_opaque_background, off",
        f"group mirror_comparison, {l_obj} {d_obj}",
        f"zoom {l_obj} or {d_obj}",
    ]
    return _write_lines(out_pml, lines)


def write_af3_correction_session(
    before_pdb: Union[str, os.PathLike],
    after_pdb: Union[str, os.PathLike],
    out_pml: Union[str, os.PathLike],
) -> Path:
    """
    Write a PyMOL script highlighting residues corrected by ChiralFold.

    Corrected residues are inferred by running ``detect_chirality_violations`` on
    ``before_pdb``. The script shows before/after structures overlaid, colors the
    uncorrected model red and corrected model green, and highlights Cbeta
    references plus corrected C-alpha centers for each changed residue.
    """
    before_report = detect_chirality_violations(str(before_pdb))
    corrected_residues = _violation_residues(before_report)
    before_obj = _object_name(before_pdb, "_before")
    after_obj = _object_name(after_pdb, "_after")
    before_sel = _selection_for_residues(before_obj, corrected_residues)
    after_sel = _selection_for_residues(after_obj, corrected_residues)
    lines = [
        "reinitialize",
        "bg_color white",
        f"load {_quote_path(before_pdb)}, {before_obj}",
        f"load {_quote_path(after_pdb)}, {after_obj}",
        f"align {after_obj}, {before_obj}",
        f"hide everything, {before_obj} or {after_obj}",
        f"show cartoon, {before_obj} or {after_obj}",
        f"show sticks, ({before_obj} or {after_obj}) and not elem H",
        f"color red, {before_obj}",
        f"color forest, {after_obj}",
        f"set cartoon_transparency, 0.45, {before_obj}",
        f"select corrected_before, {before_sel}",
        f"select corrected_after, {after_sel}",
        "show spheres, (corrected_before or corrected_after) and name CA+CB",
        "color orange, corrected_before and name CA+CB",
        "color lime, corrected_after and name CA+CB",
        'label corrected_after and name CA, "corrected %s%s" % (resn, resi)',
        "set sphere_scale, 0.35",
        "set ray_opaque_background, off",
        f"group af3_correction, {before_obj} {after_obj} corrected_before corrected_after",
        f"zoom {before_obj} or {after_obj}",
    ]
    return _write_lines(out_pml, lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Write ChiralFold PyMOL .pml scripts")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("audit", help="Color a PDB by chirality audit status")
    p.add_argument("pdb")
    p.add_argument("report_json")
    p.add_argument("-o", "--output", required=True)

    p = sub.add_parser("mirror", help="Create an L/D mirror comparison script")
    p.add_argument("l_pdb")
    p.add_argument("d_pdb")
    p.add_argument("-o", "--output", required=True)

    p = sub.add_parser("af3-correction", help="Create a before/after AF3 correction script")
    p.add_argument("before_pdb")
    p.add_argument("after_pdb")
    p.add_argument("-o", "--output", required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command == "audit":
        path = write_chirality_audit_session(args.pdb, args.report_json, args.output)
    elif args.command == "mirror":
        path = write_mirror_comparison_session(args.l_pdb, args.d_pdb, args.output)
    else:
        path = write_af3_correction_session(args.before_pdb, args.after_pdb, args.output)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
