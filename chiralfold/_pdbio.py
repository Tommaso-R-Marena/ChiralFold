"""Canonical fixed-column PDB record reader shared by the ChiralFold modules.

Before this module existed, :mod:`chiralfold.auditor`, :mod:`chiralfold.af3_correct`,
:mod:`chiralfold.interface_scorer` and :mod:`chiralfold.rotamers` each carried
their own copy of the ``ATOM``/``HETATM`` column parser. The copies had drifted
apart and disagreed on three things that change results:

* **Multi-model files.** The auditor stopped at the first ``ENDMDL``; the other
  three read every model and then silently collapsed them by atom-name key, so
  an NMR ensemble was analysed as one structure with *n* superimposed copies.
* **Alternate locations.** Three modules hard-coded the accept-list
  ``{' ', 'A', '1'}``, which drops every atom of a residue modelled only as
  altloc ``B``/``C``, and can keep a minority conformer or mix atoms from
  different alternates within one residue.
* **Insertion codes.** Two modules keyed their de-duplication on
  ``(chain, resseq, name)``, so residues ``47`` and ``47A`` overwrote each other.

Everything here is primitive-typed on purpose: each consumer still builds its own
lightweight atom object (they need different extra fields), but they all agree on
which records exist and which alternate wins.

References:
    wwPDB (2011) Protein Data Bank Contents Guide, ATOM/HETATM record format.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, Iterator, List, NamedTuple, Optional, Set, Tuple

#: Solvent residue names skipped by every ChiralFold parser.
WATER_RESNAMES: frozenset = frozenset({"HOH", "WAT", "DOD"})

#: Altloc field values meaning "no alternate location".
BLANK_ALTLOCS: frozenset = frozenset({" ", ""})


class AtomRecord(NamedTuple):
    """One parsed ``ATOM``/``HETATM`` record.

    Attributes:
        line_idx: Zero-based index of the source line (for in-place rewriting).
        record: ``"ATOM"`` or ``"HETATM"``.
        serial: Atom serial number (columns 7–11).
        name: Atom name, stripped.
        altloc: Raw alternate-location character (columns 17).
        resname: Residue name, stripped.
        chain: Chain identifier (column 22).
        resseq: Residue sequence number.
        icode: Insertion code (column 27).
        x, y, z: Orthogonal coordinates in Å.
        occupancy: Occupancy (defaults to 1.0 when absent or unparseable).
        element: Element symbol, stripped (may be empty).
        model: 1-based model number (``MODEL`` record order; 1 when absent).
    """

    line_idx: int
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    element: str
    model: int = 1

    @property
    def residue_key(self) -> Tuple[int, str, int, str]:
        """``(model, chain, resseq, icode)`` — identifies a residue."""
        return (self.model, self.chain, self.resseq, self.icode)

    @property
    def atom_key(self) -> Tuple[int, str, int, str, str]:
        """``(model, chain, resseq, icode, name)`` — identifies one atom."""
        return (self.model, self.chain, self.resseq, self.icode, self.name)


def iter_atom_records(
    lines: Iterable[str],
    *,
    first_model_only: bool = True,
    skip_water: bool = True,
) -> Iterator[AtomRecord]:
    """Yield :class:`AtomRecord` for every parseable coordinate line.

    Args:
        lines: Iterable of raw PDB lines (a file handle works).
        first_model_only: Stop at the first ``ENDMDL`` when the file is a
            multi-model (NMR / ensemble) deposition. Coordinate analysis is
            defined per model, so mixing models is never correct.
        skip_water: Drop ``HOH``/``WAT``/``DOD`` residues.

    Yields:
        Records in file order. Malformed lines and lines shorter than 54
        characters are skipped silently, matching PDB reader convention.
    """
    saw_model = False
    model = 1
    for idx, line in enumerate(lines):
        if line.startswith("MODEL"):
            if saw_model:
                if first_model_only:
                    return
                model += 1
            saw_model = True
            continue
        if line.startswith("ENDMDL"):
            if saw_model and first_model_only:
                return
            continue
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            continue

        resname = line[17:20].strip()
        if skip_water and resname in WATER_RESNAMES:
            continue

        try:
            serial = int(line[6:11])
            resseq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue

        try:
            occupancy = float(line[54:60]) if len(line) >= 60 else 1.0
        except ValueError:
            occupancy = 1.0

        yield AtomRecord(
            line_idx=idx,
            record=line[0:6].strip(),
            serial=serial,
            name=line[12:16].strip(),
            altloc=line[16],
            resname=resname,
            chain=line[21] if len(line) > 21 else " ",
            resseq=resseq,
            icode=line[26] if len(line) > 26 else " ",
            x=x, y=y, z=z,
            occupancy=occupancy,
            element=line[76:78].strip() if len(line) >= 78 else "",
            model=model,
        )


def winning_altlocs(
    records: Iterable[AtomRecord],
) -> Dict[Tuple[int, str, int, str], str]:
    """Pick one alternate-location label per residue.

    The label with the highest occupancy summed over the residue's atoms wins;
    ties break on the label itself so the choice is deterministic. Residues with
    no lettered alternates are absent from the result.

    Choosing per *residue* rather than per *atom* matters: a residue whose
    backbone came from alternate ``A`` and whose side chain came from ``B`` is
    not a conformation that exists in the model, and scoring it produces
    intra-residue geometry outliers and clashes that are pure artefacts.
    """
    weights: Dict[Tuple[int, str, int, str], Dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    for rec in records:
        if rec.altloc not in BLANK_ALTLOCS:
            weights[rec.residue_key][rec.altloc] += rec.occupancy
    return {
        rkey: min(labels.items(), key=lambda kv: (-kv[1], kv[0]))[0]
        for rkey, labels in weights.items()
    }


#: Altloc labels the pre-v3.6.0 per-atom filter accepted.
LEGACY_ALTLOC_ACCEPT: frozenset = frozenset({" ", "", "A", "1"})

#: Supported alternate-location selection policies.
ALTLOC_POLICIES: Tuple[str, ...] = ("residue", "legacy")


def select_altlocs(
    records: List[AtomRecord], policy: str = "residue"
) -> List[AtomRecord]:
    """Return *records* reduced to a single conformation.

    Args:
        records: Parsed atom records, in file order.
        policy:
            ``"residue"`` (default) — blank-altloc atoms are always kept, and
            among lettered alternates only the residue's winning label
            (see :func:`winning_altlocs`) survives. Each
            ``(model, chain, resseq, icode, name)`` appears at most once.

            ``"legacy"`` — reproduces the pre-v3.6.0 per-atom filter, which
            accepted only the labels in :data:`LEGACY_ALTLOC_ACCEPT` and then
            de-duplicated by first-seen. Provided so the effect of the policy
            change on any published number can be measured directly rather than
            argued about; see ``benchmarks/altloc_policy_sensitivity.py``.

    Raises:
        ValueError: if *policy* is not one of :data:`ALTLOC_POLICIES`.
    """
    if policy not in ALTLOC_POLICIES:
        raise ValueError(
            f"unknown altloc policy {policy!r}; expected one of {ALTLOC_POLICIES}"
        )

    winner = winning_altlocs(records) if policy == "residue" else {}
    kept: List[AtomRecord] = []
    seen: Set[Tuple[int, str, int, str, str]] = set()
    for rec in records:
        if policy == "residue":
            if rec.altloc not in BLANK_ALTLOCS:
                if winner.get(rec.residue_key) != rec.altloc:
                    continue
        elif rec.altloc not in LEGACY_ALTLOC_ACCEPT:
            continue
        key = rec.atom_key
        if key in seen:
            continue
        seen.add(key)
        kept.append(rec)
    return kept


def read_atom_records(
    path: str,
    *,
    first_model_only: bool = True,
    skip_water: bool = True,
    resolve_altlocs: bool = True,
    altloc_policy: str = "residue",
    chains: Optional[Iterable[str]] = None,
) -> List[AtomRecord]:
    """Read *path* and return its coordinate records.

    Args:
        path: PDB file path.
        first_model_only: See :func:`iter_atom_records`.
        skip_water: See :func:`iter_atom_records`.
        resolve_altlocs: Apply :func:`select_altlocs`. Pass ``False`` to keep
            every alternate (e.g. when rewriting a file verbatim).
        altloc_policy: Selection policy passed to :func:`select_altlocs`.
        chains: Optional chain-ID whitelist, applied before altloc resolution.

    Returns:
        List of :class:`AtomRecord` in file order.
    """
    with open(path) as fh:
        records = list(
            iter_atom_records(
                fh, first_model_only=first_model_only, skip_water=skip_water
            )
        )
    if chains is not None:
        allowed = set(chains)
        records = [r for r in records if r.chain in allowed]
    if resolve_altlocs:
        records = select_altlocs(records, policy=altloc_policy)
    return records


def read_lines_and_records(
    path: str,
    *,
    first_model_only: bool = True,
    skip_water: bool = True,
    resolve_altlocs: bool = True,
    altloc_policy: str = "residue",
) -> Tuple[List[str], List[AtomRecord]]:
    """Like :func:`read_atom_records` but also returns the raw lines.

    Used by writers that rewrite coordinates in place and therefore need both
    the parsed records (with ``line_idx``) and the original file content.
    """
    with open(path) as fh:
        lines = fh.readlines()
    records = list(
        iter_atom_records(
            lines, first_model_only=first_model_only, skip_water=skip_water
        )
    )
    if resolve_altlocs:
        records = select_altlocs(records, policy=altloc_policy)
    return lines, records
