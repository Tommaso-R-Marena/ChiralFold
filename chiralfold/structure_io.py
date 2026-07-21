"""Structure and sequence I/O helpers (PDB, mmCIF, FASTA).

ChiralFold's auditor / corrector / mirror pipeline reads fixed-column PDB
ATOM/HETATM records. mmCIF inputs are converted via the ``_atom_site`` loop
(no gemmi dependency). FASTA inputs are parsed for sequence + UniProt IDs
used by AlphaFold DB fetch.
"""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

# A-Z then 0-9 → 36 safe single-character PDB chain IDs.
PDB_CHAIN_CHARS = [chr(c) for c in range(ord("A"), ord("Z") + 1)] + [
    str(d) for d in range(10)
]

CONVERTER_VERSION = 2
CONVERTER_REMARK_PREFIX = "REMARK 999 CHIRALFOLD-CIF-CONVERT v"

STRUCTURE_SUFFIXES = (".pdb", ".ent", ".cif", ".mmcif")
FASTA_SUFFIXES = (".fasta", ".fa", ".faa", ".fna")
SUPPORTED_UPLOAD_SUFFIXES = STRUCTURE_SUFFIXES + FASTA_SUFFIXES

# UniProt accession (Swiss-Prot / TrEMBL style).
_UNIPROT_RE = re.compile(
    r"(?:^|[^A-Z0-9])"
    r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
    r"(?:[^A-Z0-9]|$)",
    re.IGNORECASE,
)


def mmcif_to_pdb(cif_text: str) -> Optional[str]:
    """Convert the ``_atom_site`` loop of an mmCIF document to PDB ATOM records.

    Chain IDs are remapped deterministically to single-character PDB IDs
    (A–Z then 0–9). Returns ``None`` if no atom_site loop is found or if the
    structure has more than 36 distinct chains.
    """
    lines = cif_text.splitlines()
    cols: dict = {}
    data_rows: List[List[str]] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j = i + 1
            headers: List[str] = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                headers.append(lines[j].strip())
                j += 1
            if headers and headers[0].startswith("_atom_site."):
                cols = {h.split(".", 1)[1]: k for k, h in enumerate(headers)}
                k = j
                while (
                    k < len(lines)
                    and lines[k].strip()
                    and not lines[k].strip().startswith("_")
                    and lines[k].strip() != "loop_"
                    and not lines[k].startswith("#")
                ):
                    parts = lines[k].split()
                    if len(parts) >= len(headers):
                        data_rows.append(parts)
                    k += 1
                break
            i = j
        else:
            i += 1
    if not cols or not data_rows:
        return None

    def col(parts: List[str], key: str, default: str = "") -> str:
        idx = cols.get(key)
        if idx is None or idx >= len(parts):
            return default
        val = parts[idx]
        return default if val in (".", "?") else val

    def _raw_chain(parts: List[str]) -> str:
        return col(parts, "auth_asym_id") or col(parts, "label_asym_id") or "A"

    first_model = None
    seen_order: List[str] = []
    seen_set: set = set()
    for parts in data_rows:
        model = col(parts, "pdbx_PDB_model_num", "1")
        if first_model is None:
            first_model = model
        elif model != first_model:
            break
        raw = _raw_chain(parts)
        if raw not in seen_set:
            seen_set.add(raw)
            seen_order.append(raw)
    if len(seen_order) > len(PDB_CHAIN_CHARS):
        return None
    chain_map = {mmcif_id: PDB_CHAIN_CHARS[i] for i, mmcif_id in enumerate(seen_order)}

    out_lines: List[str] = []
    serial = 0
    last_model = None
    for parts in data_rows:
        model = col(parts, "pdbx_PDB_model_num", "1")
        if last_model is not None and model != last_model:
            break
        last_model = model
        record = col(parts, "group_PDB", "ATOM")
        atom_name = col(parts, "auth_atom_id") or col(parts, "label_atom_id")
        alt = col(parts, "label_alt_id")
        resname = col(parts, "auth_comp_id") or col(parts, "label_comp_id")
        raw_chain = _raw_chain(parts)
        chain = chain_map[raw_chain]
        resseq = col(parts, "auth_seq_id") or col(parts, "label_seq_id") or "0"
        icode = col(parts, "pdbx_PDB_ins_code")
        try:
            x = float(col(parts, "Cartn_x", "0"))
            y = float(col(parts, "Cartn_y", "0"))
            z = float(col(parts, "Cartn_z", "0"))
        except ValueError:
            continue
        try:
            occ = float(col(parts, "occupancy", "1.0"))
        except ValueError:
            occ = 1.0
        try:
            bfac = float(col(parts, "B_iso_or_equiv", "0.0"))
        except ValueError:
            bfac = 0.0
        element = col(parts, "type_symbol")
        try:
            resseq_int = int(resseq)
        except ValueError:
            resseq_int = 0
        serial += 1
        if len(atom_name) >= 4:
            name_field = atom_name[:4]
        elif len(element) == 1 and not atom_name[:1].isdigit():
            name_field = " " + atom_name.ljust(3)
        else:
            name_field = atom_name.ljust(4)
        out_lines.append(
            f"{record:<6}{serial % 100000:>5} {name_field:<4}"
            f"{alt[:1]:1}{resname[:3]:>3} {chain:1}{resseq_int % 10000:>4}"
            f"{icode[:1]:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occ:6.2f}{bfac:6.2f}          {element[:2]:>2}"
        )
    if not out_lines:
        return None
    out_lines.append("END")
    header = f"{CONVERTER_REMARK_PREFIX}{CONVERTER_VERSION}"
    return header + "\n" + "\n".join(out_lines) + "\n"


def write_mmcif_as_pdb(cif_path: str, pdb_path: Optional[str] = None) -> str:
    """Convert an mmCIF file to a PDB path (temp file if *pdb_path* is None)."""
    text = Path(cif_path).read_text(encoding="utf-8", errors="replace")
    pdb_text = mmcif_to_pdb(text)
    if not pdb_text:
        raise ValueError(
            f"Could not convert mmCIF to PDB (missing _atom_site or >36 chains): "
            f"{cif_path}"
        )
    if pdb_path is None:
        fd, pdb_path = tempfile.mkstemp(prefix="chiralfold_cif_", suffix=".pdb")
        os.close(fd)
    Path(pdb_path).write_text(pdb_text)
    return pdb_path


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """Parse FASTA text into ``[(header, sequence), ...]`` (sequence upper AA)."""
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                seq = "".join(chunks).upper()
                seq = re.sub(r"[^A-Z]", "", seq)
                if seq:
                    records.append((header, seq))
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                header = "sequence"
            chunks.append(line)
    if header is not None:
        seq = "".join(chunks).upper()
        seq = re.sub(r"[^A-Z]", "", seq)
        if seq:
            records.append((header, seq))
    return records


def parse_fasta_file(path: str) -> List[Tuple[str, str]]:
    return parse_fasta(Path(path).read_text(encoding="utf-8", errors="replace"))


def uniprot_from_text(text: str) -> Optional[str]:
    """Extract the first UniProt accession from free text / FASTA header."""
    m = _UNIPROT_RE.search(text.upper())
    return m.group(1).upper() if m else None


def detect_format(path: str) -> str:
    """Return ``pdb``, ``mmcif``, ``fasta``, or ``unknown`` from suffix / peek."""
    suf = Path(path).suffix.lower()
    if suf in (".pdb", ".ent"):
        return "pdb"
    if suf in (".cif", ".mmcif"):
        return "mmcif"
    if suf in FASTA_SUFFIXES:
        return "fasta"
    # Peek for content-type when Gradio renames uploads oddly.
    try:
        head = Path(path).read_text(encoding="utf-8", errors="replace")[:400]
    except OSError:
        return "unknown"
    stripped = head.lstrip()
    if stripped.startswith(">"):
        return "fasta"
    if "data_" in stripped[:80] or "_atom_site." in head:
        return "mmcif"
    if "ATOM" in head or "HETATM" in head:
        return "pdb"
    return "unknown"


def ensure_pdb_path(path: str) -> str:
    """Return a PDB path for *path*, converting mmCIF when needed.

    FASTA is not a structure format — callers that want AFDB resolution should
    use :func:`chiralfold.fetch.resolve_to_pdb` instead.
    """
    kind = detect_format(path)
    if kind == "pdb":
        return path
    if kind == "mmcif":
        return write_mmcif_as_pdb(path)
    if kind == "fasta":
        raise ValueError(
            "FASTA has no coordinates. Provide a UniProt/AFDB ID to fetch from "
            "AlphaFold DB, or use `chiralfold predict` for sequence construction."
        )
    raise ValueError(f"Unsupported structure format: {path}")
