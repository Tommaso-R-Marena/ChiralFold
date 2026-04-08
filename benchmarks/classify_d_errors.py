#!/usr/bin/env python3
"""
classify_d_errors.py — Reproducible error taxonomy for 21 D-label/L-coordinate
mismatches in the PDB.

Downloads each raw PDB file from RCSB, extracts deposition metadata, COMPND,
SEQRES, HET/HETATM records, and computes the signed tetrahedron volume at
each Cα.  Cross-references biological context to classify each mismatch into
one of four categories:

  Stereochem  — biology requires D, coordinates show L  (most concerning)
  CCD-Code    — non-polymer ligand with wrong D-form CCD code
  Mislabel    — polymer residue labeled D in an L-protein / L-specified context
  Borderline  — inconclusive (near-zero volume, racemic conditions, etc.)

Outputs:
  results/error_taxonomy.csv          — one row per error structure (12 rows)
  results/error_taxonomy_expanded.csv — one row per individual error residue (21+ rows)
  results/error_classification.json   — summary JSON

Dependencies: numpy only (no ChiralFold code).
"""
import csv, json, os, sys, urllib.request
import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
PDB_CACHE   = RESULTS_DIR / "pdb_cache"
PDB_CACHE.mkdir(parents=True, exist_ok=True)

# Each entry defines one error-containing structure with full evidence.
ERROR_STRUCTURES = [
    {
        "pdb_id": "1ABI",
        "residue": "DPN",
        "chain": "I",
        "positions": [56],
        "error_type": "Stereochem",
        "polymer_status": "In SEQRES (polymer)",
        "record_type": "HETATM",
        "biological_context": "Hirulog-3 thrombin inhibitor peptide; D-Phe (DPN) at position 2 is by design",
        "evidence_compnd": "COMPND lists 'HIRULOG 3' peptide on chain I",
        "evidence_seqres": "DPN appears in SEQRES for chain I",
        "evidence_literature": "Hirulog peptides contain D-Phe as a key pharmacophore element",
        "evidence_internal_control": "DPN:1 vol=-2.60 (correct D); DPN:56 vol=+2.49 (L-error); same atom naming, no ALTLOCs",
        "correct_label": "DPN (D-Phe) — label is correct; coordinates are wrong",
        "conclusion": "Genuine coordinate-level stereochemistry error",
    },
    {
        "pdb_id": "1BG0",
        "residue": "DAR",
        "chain": "A",
        "positions": [403],
        "error_type": "CCD-Code",
        "polymer_status": "Not in SEQRES (non-polymer ligand)",
        "record_type": "HETATM",
        "biological_context": "Arginine kinase — substrate is L-arginine",
        "evidence_compnd": "COMPND: 'ARGININE KINASE'; title: 'TRANSITION STATE STRUCTURE OF ARGININE KINASE'",
        "evidence_seqres": "DAR does NOT appear in any SEQRES record",
        "evidence_literature": "Arginine kinase catalyses phosphoryl transfer to L-arginine (EC 2.7.3.3)",
        "evidence_internal_control": "N/A — standalone ligand, no polymer D-residues",
        "correct_label": "ARG (L-Arg) — should use L-form CCD code",
        "conclusion": "CCD code misassignment: L-arginine substrate labeled with D-Arg code (DAR)",
    },
    {
        "pdb_id": "1D7T",
        "residue": "DTY",
        "chain": "A",
        "positions": [4],
        "error_type": "Stereochem",
        "polymer_status": "In SEQRES (polymer)",
        "record_type": "HETATM",
        "biological_context": "Engineered contryphan cyclic peptide; explicitly named [D-Tyr4, Asn5, Lys7]contryphan-R",
        "evidence_compnd": "COMPND: 'YNK-CONTRYPHAN'; SEQRES: GLY CYS HYP DTY ASN PRO LYS CY3",
        "evidence_seqres": "DTY at position 4 in SEQRES for chain A",
        "evidence_literature": "Pallaghy & Norton (2000); contryphan motif CPxXPXC uses D-residue at lowercase x position",
        "evidence_internal_control": "N/A — single DTY residue in structure (NMR, B-factors all 0.00)",
        "correct_label": "DTY (D-Tyr) — label is correct; coordinates are wrong",
        "conclusion": "Genuine coordinate-level stereochemistry error in NMR structure",
    },
    {
        "pdb_id": "1HHZ",
        "residue": "DAL",
        "chain": "E",
        "positions": [1],
        "error_type": "Stereochem",
        "polymer_status": "In SEQRES (polymer)",
        "record_type": "HETATM",
        "biological_context": "Deglucobalhimycin complexed with cell wall pentapeptide (L-Ala-γ-D-Glu-L-Lys-D-Ala-D-Ala); D-Ala biologically required",
        "evidence_compnd": "COMPND: 'DEGLUCOBALHIMYCIN' + 'CELL WALL PENTAPEPTIDE' on chain E",
        "evidence_seqres": "DAL appears in SEQRES for chain E (pentapeptide)",
        "evidence_literature": "Bacterial cell wall peptidoglycan terminates in D-Ala-D-Ala; vancomycin-class antibiotics bind this motif",
        "evidence_internal_control": "0.99 Å atomic resolution — electron density unambiguous; coordinates clearly L",
        "correct_label": "DAL (D-Ala) — label is correct; coordinates are wrong",
        "conclusion": "Genuine coordinate-level stereochemistry error at atomic resolution",
    },
    {
        "pdb_id": "1KO0",
        "residue": "DLY",
        "chain": "A",
        "positions": [542],
        "error_type": "Borderline",
        "polymer_status": "Not in SEQRES (non-polymer ligand)",
        "record_type": "HETATM",
        "biological_context": "Diaminopimelate decarboxylase complexed with D,L-lysine (racemic mixture)",
        "evidence_compnd": "Title: 'CRYSTAL STRUCTURE OF A D,L-LYSINE COMPLEX OF DIAMINOPIMELATE DECARBOXYLASE'",
        "evidence_seqres": "DLY does NOT appear in SEQRES",
        "evidence_literature": "Title explicitly states D,L-lysine (racemic); could be either enantiomer in binding site",
        "evidence_internal_control": "ALTLOC B, B-factor 32.3 Å²; signed volume +0.12 (near-zero, inconclusive)",
        "correct_label": "Ambiguous — racemic conditions; could be LYS (L-Lys) from the racemic mix",
        "conclusion": "Borderline: near-zero signed volume under racemic crystallization conditions",
    },
    {
        "pdb_id": "1MCB",
        "residue": "DHI",
        "chain": "P",
        "positions": [3],
        "error_type": "Mislabel",
        "polymer_status": "In SEQRES (polymer)",
        "record_type": "HETATM",
        "biological_context": "Designed peptide ligand for immunoglobulin MCG; sequence is N-Acetyl-L-Gln-D-Phe-L-His-D-Pro-OH",
        "evidence_compnd": "COMPND explicitly states: 'PEPTIDE N-ACETYL-L-GLN-D-PHE-L-HIS-D-PRO-OH'",
        "evidence_seqres": "SEQRES: ACE GLN DPN DHI DPR — uses DHI but COMPND says L-HIS at position 3",
        "evidence_literature": "Edmundson et al. (1993) Proteins 16:246-267; peptide design with specific L/D pattern",
        "evidence_internal_control": "DPN at pos 2 (D-Phe, intended D) and DPR at pos 4 (D-Pro, intended D) are correctly specified",
        "correct_label": "HIS (L-His) — COMPND says L-HIS; DHI is wrong residue code",
        "conclusion": "Polymer residue mislabel: L-His labeled as DHI (D-His) in SEQRES despite COMPND saying L-HIS",
    },
    {
        "pdb_id": "1OF6",
        "residue": "DTY",
        "chain": "A-H",
        "positions": [1370, 1370, 1371, 1370, 1370, 1370, 1369, 1369],
        "error_type": "CCD-Code",
        "polymer_status": "Not in SEQRES (non-polymer ligand)",
        "record_type": "HETATM",
        "biological_context": "Tyrosine-regulated DAHP synthase (Aro4p) from S. cerevisiae; allosterically inhibited by L-tyrosine",
        "evidence_compnd": "Title: 'TYROSINE-REGULATED ... COMPLEXED WITH TYROSINE AND MANGANESE' (not 'D-tyrosine')",
        "evidence_seqres": "DTY does NOT appear in any SEQRES record; 0 ATOM lines, 104 HETATM lines",
        "evidence_literature": "Nature Comms Chem 2023 (doi:10.1038/s42004-023-00946-x) cites 1OF6 as Aro4 L-Tyr complex; Aro4p Ki(L-Tyr)=0.9μM",
        "evidence_internal_control": "All 8 chains: DTY has OXT (free carboxylate, standalone ligand); all volumes +2.51 to +2.67",
        "correct_label": "TYR (L-Tyr) — allosteric inhibitor is L-tyrosine, not D-tyrosine",
        "conclusion": "CCD code misassignment: L-tyrosine inhibitor labeled with D-Tyr code (DTY) in 8 copies",
    },
    {
        "pdb_id": "1P52",
        "residue": "DAR",
        "chain": "A",
        "positions": [403],
        "error_type": "CCD-Code",
        "polymer_status": "Not in SEQRES (non-polymer ligand)",
        "record_type": "HETATM",
        "biological_context": "Arginine kinase E314D mutant — same enzyme as 1BG0; substrate is L-arginine",
        "evidence_compnd": "Title: 'STRUCTURE OF ARGININE KINASE E314D MUTANT'",
        "evidence_seqres": "DAR does NOT appear in any SEQRES record",
        "evidence_literature": "Arginine kinase (EC 2.7.3.3) uses L-arginine as its physiological substrate",
        "evidence_internal_control": "N/A — standalone ligand",
        "correct_label": "ARG (L-Arg) — should use L-form CCD code",
        "conclusion": "CCD code misassignment: L-arginine substrate labeled with D-Arg code (DAR)",
    },
    {
        "pdb_id": "1UHG",
        "residue": "DSN",
        "chain": "D",
        "positions": [164],
        "error_type": "Mislabel",
        "polymer_status": "In SEQRES (polymer)",
        "record_type": "HETATM",
        "biological_context": "S-ovalbumin; standard L-protein from chicken egg white — no biological reason for D-Ser",
        "evidence_compnd": "Title: 'CRYSTAL STRUCTURE OF S-OVALBUMIN AT 1.9 ANGSTROM RESOLUTION'",
        "evidence_seqres": "DSN appears in SEQRES — but ovalbumin is an L-protein (UniProt P01012)",
        "evidence_literature": "Ovalbumin contains only standard L-amino acids; no known post-translational D-amino acid modifications",
        "evidence_internal_control": "All other residues in ovalbumin are standard L-forms",
        "correct_label": "SER (L-Ser) — standard L-protein; DSN is wrong residue code",
        "conclusion": "Polymer residue mislabel: L-Ser at position 164 labeled as DSN (D-Ser) in an L-protein",
    },
    {
        "pdb_id": "1XT7",
        "residue": "DSG",
        "chain": "A",
        "positions": [3],
        "error_type": "Stereochem",
        "polymer_status": "In SEQRES (polymer)",
        "record_type": "HETATM",
        "biological_context": "Daptomycin (cyclic lipopeptide antibiotic); contains D-Asn (DSG) by biosynthetic requirement",
        "evidence_compnd": "COMPND: 'DAPTOMYCIN'; NMR structure",
        "evidence_seqres": "DSG appears in SEQRES for chain A",
        "evidence_literature": "Daptomycin structure: D-Asn at position 2 is required for antimicrobial activity and ring closure",
        "evidence_internal_control": "NMR structure (B-factors 0.00); other D-residues in daptomycin should be checked",
        "correct_label": "DSG (D-Asn) — label is correct; coordinates are wrong",
        "conclusion": "Genuine coordinate-level stereochemistry error in daptomycin NMR structure",
    },
    {
        "pdb_id": "2AOU",
        "residue": "DCY",
        "chain": "A",
        "positions": [248],
        "error_type": "Mislabel",
        "polymer_status": "In SEQRES (polymer)",
        "record_type": "HETATM",
        "biological_context": "Histamine N-methyltransferase complexed with amodiaquine; standard L-protein (human enzyme)",
        "evidence_compnd": "Title: 'HISTAMINE METHYLTRANSFERASE COMPLEXED WITH THE ANTIMALARIAL DRUG AMODIAQUINE'",
        "evidence_seqres": "DCY appears in SEQRES — but HNMT is an L-protein (UniProt P50135)",
        "evidence_literature": "Human histamine N-methyltransferase contains only standard L-amino acids",
        "evidence_internal_control": "All other residues are standard L-forms; no post-translational D-modifications known",
        "correct_label": "CYS (L-Cys) — standard L-protein; DCY is wrong residue code",
        "conclusion": "Polymer residue mislabel: L-Cys at position 248 labeled as DCY (D-Cys) in an L-protein",
    },
    {
        "pdb_id": "2ATS",
        "residue": "DLY",
        "chain": "A",
        "positions": [3001, 3002, 3003],
        "error_type": "CCD-Code",
        "polymer_status": "Not in SEQRES (non-polymer ligand)",
        "record_type": "HETATM",
        "biological_context": "Dihydrodipicolinate synthase co-crystallised with (S)-lysine [= L-lysine]",
        "evidence_compnd": "Title: 'DIHYDRODIPICOLINATE SYNTHASE CO-CRYSTALLISED WITH (S)-LYSINE'",
        "evidence_seqres": "DLY does NOT appear in any SEQRES record",
        "evidence_literature": "(S)-lysine is the IUPAC designation for L-lysine; enzyme is feedback-inhibited by L-lysine",
        "evidence_internal_control": "3 copies (one per chain in trimer); all volumes +2.56 to +2.59",
        "correct_label": "LYS (L-Lys) — title says '(S)-lysine'; DLY is wrong CCD code",
        "conclusion": "CCD code misassignment: L-lysine labeled with D-Lys code (DLY) in 3 copies",
    },
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def download_pdb(pdb_id: str) -> Path:
    """Download PDB file from RCSB if not cached."""
    path = PDB_CACHE / f"{pdb_id}.pdb"
    if not path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"  Downloading {url} ...")
        urllib.request.urlretrieve(url, path)
    return path


def parse_pdb_metadata(path: Path) -> dict:
    """Extract title, deposition date, resolution, REMARK 2 from PDB file."""
    title_parts, dep_date, resolution = [], "", ""
    with open(path) as f:
        for line in f:
            if line.startswith("TITLE"):
                title_parts.append(line[10:].strip())
            if line.startswith("REVDAT   1"):
                dep_date = line[13:22].strip()
            if line.startswith("HEADER"):
                dep_date = line[50:59].strip()
            if line.startswith("REMARK   2 RESOLUTION."):
                res_str = line[22:].strip().replace("ANGSTROMS.", "").strip()
                if res_str and res_str != "NOT APPLICABLE.":
                    resolution = res_str + " Å"
                else:
                    resolution = "NMR"
    return {
        "title": " ".join(title_parts),
        "deposition_date": dep_date,
        "resolution": resolution,
    }


def check_seqres(path: Path, resname: str) -> bool:
    """Check if residue name appears in any SEQRES record."""
    with open(path) as f:
        for line in f:
            if line.startswith("SEQRES") and resname in line:
                return True
    return False


def check_compnd(path: Path) -> str:
    """Extract COMPND records."""
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith("COMPND"):
                parts.append(line[10:].strip())
    return " ".join(parts)


def count_atom_hetatm(path: Path, resname: str) -> tuple:
    """Count ATOM vs HETATM lines for a given residue name."""
    n_atom, n_hetatm = 0, 0
    with open(path) as f:
        for line in f:
            if line[17:20].strip() == resname:
                if line.startswith("ATOM"):
                    n_atom += 1
                elif line.startswith("HETATM"):
                    n_hetatm += 1
    return n_atom, n_hetatm


def signed_volume(path: Path, resname: str, chain: str, resseq: int) -> float:
    """Compute signed tetrahedron volume at Cα for a specific residue."""
    atoms = {}
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[17:20].strip() != resname:
                continue
            if line[21] != chain:
                continue
            if int(line[22:26]) != resseq:
                continue
            aname = line[12:16].strip()
            if aname in ("N", "CA", "C", "CB"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms[aname] = np.array([x, y, z])
    if not all(a in atoms for a in ("N", "CA", "C", "CB")):
        return float("nan")
    v1 = atoms["N"] - atoms["CA"]
    v2 = atoms["C"] - atoms["CA"]
    v3 = atoms["CB"] - atoms["CA"]
    return float(np.dot(v3, np.cross(v1, v2)))


def get_bfactor_altloc(path: Path, resname: str, chain: str, resseq: int) -> tuple:
    """Get B-factor of CA and ALTLOC status."""
    bfac, altloc_chars = None, set()
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[17:20].strip() != resname or line[21] != chain:
                continue
            if int(line[22:26]) != resseq:
                continue
            alt = line[16].strip()
            if alt:
                altloc_chars.add(alt)
            aname = line[12:16].strip()
            if aname == "CA":
                bfac = float(line[60:66])
    bfac_str = f"{bfac:.1f}" if bfac is not None else "N/A"
    altloc_str = ",".join(sorted(altloc_chars)) if altloc_chars else "No"
    return bfac_str, altloc_str


def has_oxt(path: Path, resname: str, chain: str, resseq: int) -> bool:
    """Check if residue has OXT atom (free carboxylate = standalone ligand)."""
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[17:20].strip() != resname or line[21] != chain:
                continue
            if int(line[22:26]) != resseq:
                continue
            if line[12:16].strip() == "OXT":
                return True
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 80)
    print("D-Amino Acid Error Classification — Reproducible Dataset Builder")
    print("=" * 80)

    # Per-structure rows
    structure_rows = []
    # Per-residue rows (expanded)
    residue_rows = []

    for entry in ERROR_STRUCTURES:
        pdb_id = entry["pdb_id"]
        resname = entry["residue"]
        print(f"\n--- {pdb_id} ({resname}) ---")

        pdb_path = download_pdb(pdb_id)
        meta = parse_pdb_metadata(pdb_path)
        in_seqres = check_seqres(pdb_path, resname)
        compnd = check_compnd(pdb_path)
        n_atom, n_hetatm = count_atom_hetatm(pdb_path, resname)

        # Determine chains from entry
        chain_str = entry["chain"]
        if "-" in chain_str and len(chain_str) == 3:
            chains = [chr(c) for c in range(ord(chain_str[0]), ord(chain_str[2]) + 1)]
        else:
            chains = [chain_str]

        # Compute volumes for each position/chain
        volumes = []
        bfactors = []
        altlocs = []
        oxt_flags = []

        positions = entry["positions"]
        if len(positions) == 1 and len(chains) == 1:
            # Single residue
            vol = signed_volume(pdb_path, resname, chains[0], positions[0])
            bfac, altloc = get_bfactor_altloc(pdb_path, resname, chains[0], positions[0])
            oxt = has_oxt(pdb_path, resname, chains[0], positions[0])
            volumes.append(vol)
            bfactors.append(bfac)
            altlocs.append(altloc)
            oxt_flags.append(oxt)

            residue_rows.append({
                "PDB_ID": pdb_id,
                "Residue_Code": resname,
                "Chain": chains[0],
                "Position": positions[0],
                "Signed_Volume": f"{vol:+.2f}",
                "B_Factor_CA": bfac,
                "ALTLOC": altloc,
                "Has_OXT": oxt,
                "In_SEQRES": in_seqres,
                "ATOM_Lines": n_atom,
                "HETATM_Lines": n_hetatm,
                "Error_Type": entry["error_type"],
                "Correct_Label": entry["correct_label"],
            })
        elif len(positions) > 1 and len(chains) == 1:
            # Multiple positions in same chain (e.g. 2ATS)
            for pos in positions:
                vol = signed_volume(pdb_path, resname, chains[0], pos)
                bfac, altloc = get_bfactor_altloc(pdb_path, resname, chains[0], pos)
                oxt = has_oxt(pdb_path, resname, chains[0], pos)
                volumes.append(vol)
                bfactors.append(bfac)
                altlocs.append(altloc)
                oxt_flags.append(oxt)

                residue_rows.append({
                    "PDB_ID": pdb_id,
                    "Residue_Code": resname,
                    "Chain": chains[0],
                    "Position": pos,
                    "Signed_Volume": f"{vol:+.2f}",
                    "B_Factor_CA": bfac,
                    "ALTLOC": altloc,
                    "Has_OXT": oxt,
                    "In_SEQRES": in_seqres,
                    "ATOM_Lines": n_atom,
                    "HETATM_Lines": n_hetatm,
                    "Error_Type": entry["error_type"],
                    "Correct_Label": entry["correct_label"],
                })
        else:
            # Multiple chains (e.g. 1OF6: 8 chains with different positions)
            for ch_idx, ch in enumerate(chains):
                pos = positions[ch_idx] if ch_idx < len(positions) else positions[0]
                vol = signed_volume(pdb_path, resname, ch, pos)
                bfac, altloc = get_bfactor_altloc(pdb_path, resname, ch, pos)
                oxt = has_oxt(pdb_path, resname, ch, pos)
                volumes.append(vol)
                bfactors.append(bfac)
                altlocs.append(altloc)
                oxt_flags.append(oxt)

                residue_rows.append({
                    "PDB_ID": pdb_id,
                    "Residue_Code": resname,
                    "Chain": ch,
                    "Position": pos,
                    "Signed_Volume": f"{vol:+.2f}",
                    "B_Factor_CA": bfac,
                    "ALTLOC": altloc,
                    "Has_OXT": oxt,
                    "In_SEQRES": in_seqres,
                    "ATOM_Lines": n_atom,
                    "HETATM_Lines": n_hetatm,
                    "Error_Type": entry["error_type"],
                    "Correct_Label": entry["correct_label"],
                })

        vol_range = (
            f"{min(volumes):+.2f}"
            if len(volumes) == 1
            else f"{min(volumes):+.2f} to {max(volumes):+.2f}"
        )
        n_errors = len(volumes)

        print(f"  Title:      {meta['title'][:90]}")
        print(f"  Resolution: {meta['resolution']}")
        print(f"  Dep. date:  {meta['deposition_date']}")
        print(f"  In SEQRES:  {in_seqres}  |  ATOM: {n_atom}  HETATM: {n_hetatm}")
        print(f"  Volumes:    {vol_range}  ({n_errors} residues)")
        print(f"  Has OXT:    {any(oxt_flags)}  (standalone ligand indicator)")
        print(f"  Error type: {entry['error_type']}")
        print(f"  Conclusion: {entry['conclusion']}")

        structure_rows.append({
            "PDB_ID": pdb_id,
            "Residue_Code": resname,
            "Chain": entry["chain"],
            "Positions": ";".join(str(p) for p in sorted(set(entry["positions"]))),
            "N_Errors": n_errors,
            "Signed_Volume_Range": vol_range,
            "Resolution": meta["resolution"],
            "Deposition_Date": meta["deposition_date"],
            "Title": meta["title"],
            "In_SEQRES": in_seqres,
            "ATOM_Lines": n_atom,
            "HETATM_Lines": n_hetatm,
            "Polymer_Status": entry["polymer_status"],
            "Error_Type": entry["error_type"],
            "Biological_Context": entry["biological_context"],
            "Evidence_COMPND": entry["evidence_compnd"],
            "Evidence_SEQRES": entry["evidence_seqres"],
            "Evidence_Literature": entry["evidence_literature"],
            "Evidence_Internal_Control": entry["evidence_internal_control"],
            "Correct_Label": entry["correct_label"],
            "Conclusion": entry["conclusion"],
        })

    # -----------------------------------------------------------------------
    # Write CSV outputs
    # -----------------------------------------------------------------------
    # Per-structure taxonomy
    csv_struct = RESULTS_DIR / "error_taxonomy.csv"
    struct_fields = list(structure_rows[0].keys())
    with open(csv_struct, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=struct_fields)
        writer.writeheader()
        writer.writerows(structure_rows)
    print(f"\n\nWrote {csv_struct}  ({len(structure_rows)} structures)")

    # Per-residue expanded
    csv_expanded = RESULTS_DIR / "error_taxonomy_expanded.csv"
    exp_fields = list(residue_rows[0].keys())
    with open(csv_expanded, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=exp_fields)
        writer.writeheader()
        writer.writerows(residue_rows)
    print(f"Wrote {csv_expanded}  ({len(residue_rows)} residues)")

    # Summary JSON
    type_counts = {}
    for row in structure_rows:
        t = row["Error_Type"]
        if t not in type_counts:
            type_counts[t] = {"structures": 0, "errors": 0, "pdbs": []}
        type_counts[t]["structures"] += 1
        type_counts[t]["errors"] += row["N_Errors"]
        type_counts[t]["pdbs"].append(row["PDB_ID"])

    summary = {
        "total_mismatches": sum(tc["errors"] for tc in type_counts.values()),
        "total_structures": sum(tc["structures"] for tc in type_counts.values()),
        "classification": type_counts,
        "method": "Signed tetrahedron volume at Cα: V = (CB-CA) · ((N-CA) × (C-CA))",
        "convention": "V > 0 → L-stereochemistry; V < 0 → D-stereochemistry",
        "dependencies": "numpy only; no ChiralFold code",
        "data_source": "Raw PDB files from https://files.rcsb.org/download/{PDB_ID}.pdb",
        "classification_evidence": "COMPND records, SEQRES, HET/HETATM, deposition titles, primary literature",
    }

    json_path = RESULTS_DIR / "error_classification.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote {json_path}")

    # -----------------------------------------------------------------------
    # Print summary
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for t in ["Stereochem", "CCD-Code", "Mislabel", "Borderline"]:
        if t in type_counts:
            tc = type_counts[t]
            print(f"  {t:12s}: {tc['errors']:2d} errors in {tc['structures']} structures  {tc['pdbs']}")
    print(f"  {'TOTAL':12s}: {summary['total_mismatches']} mismatches in {summary['total_structures']} structures")
    print("=" * 80)


if __name__ == "__main__":
    main()
