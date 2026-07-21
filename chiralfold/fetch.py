"""Download structures from RCSB PDB and AlphaFold DB (when online).

Accepts:
  - 4-character PDB IDs → RCSB ``.pdb`` (falls back to ``.cif``)
  - UniProt accessions (e.g. ``P04637``) → AlphaFold DB model
  - AFDB entry tags (e.g. ``AF-P04637-F1``) → AlphaFold DB model
  - Local PDB / mmCIF / FASTA paths (FASTA → UniProt in header → AFDB)

No heavy model inference — download only.
"""

from __future__ import annotations

import os
import re
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

from .structure_io import (
    detect_format,
    ensure_pdb_path,
    parse_fasta_file,
    uniprot_from_text,
    write_mmcif_as_pdb,
)

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
AFDB_FILE_URL = "https://alphafold.ebi.ac.uk/files/{filename}"
AFDB_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"

_USER_AGENT = "ChiralFold/3.5 (+https://github.com/Tommaso-R-Marena/ChiralFold)"
_DEFAULT_TIMEOUT = 45

_PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_AFDB_RE = re.compile(
    r"^AF-([A-Z0-9]+)-F(\d+)(?:-model_v\d+)?$",
    re.IGNORECASE,
)
_UNIPROT_ONLY_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$",
    re.IGNORECASE,
)


@dataclass
class FetchResult:
    """Result of a structure fetch / local resolve."""

    path: str
    source: str  # rcsb | afdb | local | local-mmcif | fasta-afdb
    query: str
    label: str
    bytes_downloaded: int = 0


class FetchError(ValueError):
    """User-facing fetch / resolve failure."""


def _http_get(url: str, timeout: float = _DEFAULT_TIMEOUT) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def classify_query(query: str) -> Tuple[str, str]:
    """Return ``(kind, normalized)`` where kind is pdb|uniprot|afdb|path|unknown."""
    q = (query or "").strip()
    if not q:
        return "unknown", ""
    if os.path.isfile(q):
        return "path", q
    # Strip common wrappers
    q_clean = q.strip().strip("\"'")
    af = _AFDB_RE.match(q_clean)
    if af:
        return "afdb", f"AF-{af.group(1).upper()}-F{af.group(2)}"
    if _PDB_ID_RE.match(q_clean):
        return "pdb", q_clean.upper()
    if _UNIPROT_ONLY_RE.match(q_clean):
        return "uniprot", q_clean.upper()
    # AFDB full filename without extension
    if q_clean.upper().startswith("AF-") and "-F" in q_clean.upper():
        base = q_clean
        for suf in (".pdb", ".cif", ".PDB", ".CIF"):
            if base.endswith(suf):
                base = base[: -len(suf)]
        af2 = _AFDB_RE.match(base)
        if af2:
            return "afdb", f"AF-{af2.group(1).upper()}-F{af2.group(2)}"
    up = uniprot_from_text(q_clean)
    if up and len(q_clean) <= 40:
        return "uniprot", up
    return "unknown", q_clean


def _write_temp(data: bytes, suffix: str, stem: str) -> str:
    out_dir = Path(tempfile.mkdtemp(prefix="chiralfold_fetch_"))
    path = out_dir / f"{stem}{suffix}"
    path.write_bytes(data)
    return str(path)


def _looks_like_pdb(data: bytes) -> bool:
    return b"ATOM" in data or b"HETATM" in data


def _looks_like_cif(data: bytes) -> bool:
    head = data[:800].decode("utf-8", errors="replace")
    return "data_" in head or "_atom_site." in head


def fetch_rcsb(
    pdb_id: str,
    *,
    timeout: float = _DEFAULT_TIMEOUT,
    max_bytes: Optional[int] = None,
) -> FetchResult:
    """Download a structure from RCSB (PDB, then mmCIF fallback)."""
    pid = pdb_id.upper().strip()
    if not _PDB_ID_RE.match(pid):
        raise FetchError(f"Invalid PDB ID: {pdb_id!r} (expected 4 characters).")

    pdb_err: Optional[str] = None
    try:
        data = _http_get(RCSB_PDB_URL.format(pdb_id=pid), timeout=timeout)
        if max_bytes is not None and len(data) > max_bytes:
            raise FetchError(f"{pid} exceeds size limit ({len(data)} bytes).")
        if _looks_like_pdb(data):
            path = _write_temp(data, ".pdb", pid)
            return FetchResult(
                path=path,
                source="rcsb",
                query=pid,
                label=pid,
                bytes_downloaded=len(data),
            )
        pdb_err = "downloaded .pdb lacked ATOM/HETATM records"
    except urllib.error.HTTPError as exc:
        pdb_err = f"HTTP {exc.code}"
    except FetchError:
        raise
    except Exception as exc:  # noqa: BLE001
        pdb_err = str(exc)

    try:
        data = _http_get(RCSB_CIF_URL.format(pdb_id=pid), timeout=timeout)
    except Exception as exc:  # noqa: BLE001
        raise FetchError(
            f"RCSB fetch failed for {pid} (PDB: {pdb_err}; CIF: {exc})."
        ) from exc
    if max_bytes is not None and len(data) > max_bytes:
        raise FetchError(f"{pid} mmCIF exceeds size limit ({len(data)} bytes).")
    if not _looks_like_cif(data):
        raise FetchError(f"RCSB mmCIF for {pid} does not look like mmCIF.")
    cif_path = _write_temp(data, ".cif", pid)
    try:
        pdb_path = write_mmcif_as_pdb(cif_path)
    except ValueError as exc:
        raise FetchError(str(exc)) from exc
    return FetchResult(
        path=pdb_path,
        source="rcsb",
        query=pid,
        label=f"{pid} (mmCIF→PDB)",
        bytes_downloaded=len(data),
    )


def _afdb_candidates(uniprot: str, fragment: int = 1) -> list:
    uid = uniprot.upper()
    frag = f"F{fragment}"
    names = []
    for ver in ("v4", "v3", "v2"):
        names.append(f"AF-{uid}-{frag}-model_{ver}.pdb")
    for ver in ("v4", "v3"):
        names.append(f"AF-{uid}-{frag}-model_{ver}.cif")
    return names


def fetch_afdb(
    uniprot_or_afdb: str,
    *,
    fragment: int = 1,
    timeout: float = _DEFAULT_TIMEOUT,
    max_bytes: Optional[int] = None,
) -> FetchResult:
    """Download an AlphaFold DB model by UniProt accession or AF-…-Fn tag."""
    kind, norm = classify_query(uniprot_or_afdb)
    if kind == "afdb":
        m = _AFDB_RE.match(norm)
        assert m
        uid, fragment = m.group(1).upper(), int(m.group(2))
    elif kind == "uniprot":
        uid = norm
    else:
        uid = uniprot_from_text(uniprot_or_afdb)
        if not uid:
            raise FetchError(
                f"Not a UniProt / AFDB id: {uniprot_or_afdb!r}."
            )

    errors = []
    # Prefer API pdbUrl when available.
    try:
        import json

        raw = _http_get(AFDB_API_URL.format(uniprot=uid), timeout=timeout)
        payload = json.loads(raw.decode("utf-8"))
        entries = payload if isinstance(payload, list) else [payload]
        for entry in entries:
            for key in ("pdbUrl", "cifUrl"):
                url = entry.get(key) if isinstance(entry, dict) else None
                if not url:
                    continue
                try:
                    data = _http_get(url, timeout=timeout)
                    if max_bytes is not None and len(data) > max_bytes:
                        raise FetchError(
                            f"AFDB model for {uid} exceeds size limit."
                        )
                    return _materialize_afdb(data, uid, fragment, len(data))
                except FetchError:
                    raise
                except Exception as exc:  # noqa: BLE001
                    errors.append(f"{key}: {exc}")
    except FetchError:
        raise
    except Exception as exc:  # noqa: BLE001
        errors.append(f"api: {exc}")

    for name in _afdb_candidates(uid, fragment):
        url = AFDB_FILE_URL.format(filename=name)
        try:
            data = _http_get(url, timeout=timeout)
        except urllib.error.HTTPError as exc:
            errors.append(f"{name}: HTTP {exc.code}")
            continue
        except Exception as exc:  # noqa: BLE001
            errors.append(f"{name}: {exc}")
            continue
        if max_bytes is not None and len(data) > max_bytes:
            raise FetchError(f"AFDB model for {uid} exceeds size limit.")
        return _materialize_afdb(data, uid, fragment, len(data), name)

    detail = "; ".join(errors[:4]) if errors else "not found"
    raise FetchError(f"AlphaFold DB fetch failed for {uid}: {detail}.")


def _materialize_afdb(
    data: bytes,
    uid: str,
    fragment: int,
    nbytes: int,
    name: str = "",
) -> FetchResult:
    stem = f"AF-{uid}-F{fragment}"
    as_cif = name.endswith(".cif") or (
        _looks_like_cif(data) and b"data_" in data[:240]
    )
    if as_cif:
        cif_path = _write_temp(data, ".cif", stem)
        pdb_path = write_mmcif_as_pdb(cif_path)
        return FetchResult(
            path=pdb_path,
            source="afdb",
            query=uid,
            label=stem,
            bytes_downloaded=nbytes,
        )
    if _looks_like_pdb(data):
        path = _write_temp(data, ".pdb", stem)
        return FetchResult(
            path=path,
            source="afdb",
            query=uid,
            label=stem,
            bytes_downloaded=nbytes,
        )
    raise FetchError(f"AFDB download for {uid} was not PDB or mmCIF.")


def fetch_structure(
    query: str,
    *,
    timeout: float = _DEFAULT_TIMEOUT,
    max_bytes: Optional[int] = None,
) -> FetchResult:
    """Fetch or resolve *query* to a local PDB path."""
    kind, norm = classify_query(query)
    if kind == "path":
        return resolve_local(norm, max_bytes=max_bytes, timeout=timeout)
    if kind == "pdb":
        return fetch_rcsb(norm, timeout=timeout, max_bytes=max_bytes)
    if kind in ("uniprot", "afdb"):
        return fetch_afdb(norm, timeout=timeout, max_bytes=max_bytes)
    raise FetchError(
        "Unrecognized query. Use a PDB ID (1UBQ), UniProt accession (P04637), "
        f"AFDB id (AF-P04637-F1), or a local file path. Got: {query!r}."
    )


def resolve_local(
    path: str,
    *,
    max_bytes: Optional[int] = None,
    timeout: float = _DEFAULT_TIMEOUT,
) -> FetchResult:
    """Resolve a local PDB / mmCIF / FASTA path to a PDB file."""
    p = Path(path)
    if not p.is_file():
        raise FetchError(f"File not found: {path}")
    size = p.stat().st_size
    if max_bytes is not None and size > max_bytes:
        raise FetchError(f"File too large ({size} bytes).")
    if size == 0:
        raise FetchError("File is empty.")

    kind = detect_format(path)
    if kind == "pdb":
        return FetchResult(path=str(p), source="local", query=path, label=p.name)
    if kind == "mmcif":
        pdb_path = write_mmcif_as_pdb(str(p))
        return FetchResult(
            path=pdb_path,
            source="local-mmcif",
            query=path,
            label=f"{p.name} (mmCIF→PDB)",
        )
    if kind == "fasta":
        records = parse_fasta_file(str(p))
        if not records:
            raise FetchError(f"No sequences found in FASTA: {path}")
        header, _seq = records[0]
        uid = uniprot_from_text(header) or uniprot_from_text(p.stem)
        if not uid:
            raise FetchError(
                "FASTA has no UniProt accession in the header. "
                "Add one (e.g. >sp|P04637|P53_HUMAN) or fetch by UniProt/AFDB id."
            )
        result = fetch_afdb(uid, timeout=timeout, max_bytes=max_bytes)
        result.source = "fasta-afdb"
        result.query = path
        result.label = f"{uid} via FASTA header"
        return result
    # Last resort: try ensure_pdb_path
    try:
        pdb_path = ensure_pdb_path(str(p))
    except ValueError as exc:
        raise FetchError(str(exc)) from exc
    return FetchResult(path=pdb_path, source="local", query=path, label=p.name)


def resolve_to_pdb(
    query_or_path: str,
    *,
    timeout: float = _DEFAULT_TIMEOUT,
    max_bytes: Optional[int] = None,
) -> str:
    """Convenience: return a PDB path for a query or local file."""
    return fetch_structure(
        query_or_path, timeout=timeout, max_bytes=max_bytes
    ).path
