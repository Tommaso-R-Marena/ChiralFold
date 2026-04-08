"""
ChiralFold — Command-Line Interface
=====================================

Usage:
    chiralfold predict SEQUENCE [--chirality PATTERN]
    chiralfold validate SEQUENCE --chirality PATTERN
    chiralfold audit PDB_FILE [--chain CHAIN] [--json]
    chiralfold mirror PDB_FILE [--output OUTPUT] [--chains A,B]
    chiralfold mirror-id PDB_ID [--output OUTPUT] [--chains A,B]
    chiralfold benchmark [--output DIR]
"""

import argparse
import json
import sys

from . import __version__
from .model import ChiralFold, mixed_peptide_smiles
from .validator import validate_diastereomer


def main():
    parser = argparse.ArgumentParser(
        prog='chiralfold',
        description=(
            'ChiralFold: Stereochemistry toolkit for proteins and peptides — '
            'chirality-preserving structure generation, PDB auditing, '
            'and mirror-image transformation'
        ),
    )
    parser.add_argument(
        '--version', action='version', version=f'chiralfold {__version__}'
    )

    sub = parser.add_subparsers(dest='command', help='Available commands')

    # ── predict ───────────────────────────────────────────────────────────
    p = sub.add_parser('predict', help='Predict peptide structure')
    p.add_argument('sequence', help='Amino acid sequence (one-letter codes)')
    p.add_argument('--chirality', '-c', default=None,
                   help='Per-residue chirality (D/L string). Default: all D.')
    p.add_argument('--json', action='store_true', help='JSON output')

    # ── validate ──────────────────────────────────────────────────────────
    p = sub.add_parser('validate', help='Validate chirality of a peptide')
    p.add_argument('sequence', help='Amino acid sequence')
    p.add_argument('--chirality', '-c', required=True,
                   help='Per-residue chirality (D/L string).')
    p.add_argument('--json', action='store_true', help='JSON output')

    # ── audit (NEW in v3) ─────────────────────────────────────────────────
    p = sub.add_parser(
        'audit',
        help='Audit any PDB structure for stereochemical quality'
    )
    p.add_argument('pdb_file', help='Path to PDB file')
    p.add_argument('--chain', default=None,
                   help='Specific chain to audit (default: all)')
    p.add_argument('--json', action='store_true', help='JSON output')

    # ── mirror (NEW in v3) ────────────────────────────────────────────────
    p = sub.add_parser(
        'mirror',
        help='Mirror a PDB structure (L↔D enantiomer transformation)'
    )
    p.add_argument('pdb_file', help='Path to input PDB file')
    p.add_argument('--output', '-o', default=None,
                   help='Output PDB path (default: <input>_mirror.pdb)')
    p.add_argument('--chains', default=None,
                   help='Comma-separated chain IDs to mirror (default: all)')

    # ── mirror-id (NEW in v3) ─────────────────────────────────────────────
    p = sub.add_parser(
        'mirror-id',
        help='Download PDB by ID and mirror (L↔D transformation)'
    )
    p.add_argument('pdb_id', help='4-character PDB ID (e.g. 1SHG)')
    p.add_argument('--output', '-o', default=None,
                   help='Output PDB path')
    p.add_argument('--chains', default=None,
                   help='Comma-separated chain IDs to mirror')

    # ── benchmark ─────────────────────────────────────────────────────────
    p = sub.add_parser('benchmark', help='Run full benchmark suite')
    p.add_argument('--output', '-o', default='.', help='Output directory')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    handlers = {
        'predict': _cmd_predict,
        'validate': _cmd_validate,
        'audit': _cmd_audit,
        'mirror': _cmd_mirror,
        'mirror-id': _cmd_mirror_id,
        'benchmark': _cmd_benchmark,
    }
    handlers[args.command](args)


def _cmd_predict(args):
    seq = args.sequence.upper()
    chir = (args.chirality or 'D' * len(seq)).upper()
    model = ChiralFold(n_conformers=10)
    result = model.predict(seq, chirality_pattern=chir)

    if 'error' in result:
        print(f"Error: {result['error']}", file=sys.stderr)
        sys.exit(1)

    if args.json:
        out = {k: v for k, v in result.items() if k != 'mol'}
        print(json.dumps(out, indent=2, default=str))
    else:
        print(f"Sequence:          {result['sequence']}")
        print(f"Chirality pattern: {result['chirality_pattern']}")
        print(f"D-residues:        {result['n_d_residues']}")
        print(f"L-residues:        {result['n_l_residues']}")
        print(f"SMILES:            {result['smiles']}")
        print(f"Violations:        {result['chirality_violations']} "
              f"({result['violation_rate']:.1%})")
        if result.get('conformers'):
            print(f"3D conformers:     {len(result['conformers'])}")


def _cmd_validate(args):
    seq = args.sequence.upper()
    chir = args.chirality.upper()
    report = validate_diastereomer(seq, chir)

    if args.json:
        print(json.dumps(report, indent=2, default=str))
    else:
        status = 'PASS' if report['valid'] else 'FAIL'
        print(f"[{status}] {seq} ({chir})")
        print(f"  Residues: {report['n_residues']} "
              f"(D={report['n_d']}, L={report['n_l']}, G={report['n_glycine']})")
        print(f"  SMILES violations: {report['smiles_violations']}")
        if report['geom_checked'] > 0:
            print(f"  3D geometry: {report['geom_correct']}/{report['geom_checked']} "
                  f"correct, {report['geom_planar']} planar")


def _cmd_audit(args):
    from .auditor import audit_pdb, format_report

    report = audit_pdb(args.pdb_file)

    if args.json:
        # Remove non-serializable items
        safe = {k: v for k, v in report.items()
                if k not in ('_atoms',)}
        print(json.dumps(safe, indent=2, default=str))
    else:
        format_report(report)


def _cmd_mirror(args):
    from .pdb_pipeline import mirror_pdb

    output = args.output
    if output is None:
        base = args.pdb_file.rsplit('.', 1)[0]
        output = f"{base}_mirror.pdb"

    chains = args.chains.split(',') if args.chains else None
    result = mirror_pdb(args.pdb_file, output, chains=chains)

    print(f"Mirrored {result['n_atoms']} atoms, "
          f"{result['n_residues']} residues")
    print(f"Output: {result['output_path']}")
    print(f"Renamed: {result['stats']['residues_renamed']} residue records")


def _cmd_mirror_id(args):
    from .pdb_pipeline import fetch_and_mirror

    output = args.output or f"{args.pdb_id.upper()}_mirror.pdb"
    chains = args.chains.split(',') if args.chains else None
    result = fetch_and_mirror(args.pdb_id, output, chains=chains)

    print(f"Downloaded {args.pdb_id.upper()} from RCSB")
    print(f"Mirrored {result['n_atoms']} atoms, "
          f"{result['n_residues']} residues")
    print(f"Output: {result['output_path']}")


def _cmd_benchmark(args):
    from .data.test_sequences import PURE_D_SEQS, DIASTEREOMER_SEQS
    from .validator import validate_smiles_chirality
    from .model import d_peptide_smiles
    from rdkit import Chem

    print("Running ChiralFold benchmark suite...")
    total_chiral = total_viol = 0

    print("\nPure D-peptides:")
    for sid, seq in PURE_D_SEQS.items():
        smi = d_peptide_smiles(seq)
        mol = Chem.MolFromSmiles(smi)
        sv = validate_smiles_chirality(mol, seq, 'D' * len(seq))
        total_chiral += sv['n_chiral']
        total_viol += sv['violations']
        print(f"  {sid:<8} {seq:<22} chiral={sv['n_chiral']:>2} viol={sv['violations']}")

    print("\nDiastereomers:")
    for sid, data in DIASTEREOMER_SEQS.items():
        smi = mixed_peptide_smiles(data['seq'], data['chirality'])
        mol = Chem.MolFromSmiles(smi)
        sv = validate_smiles_chirality(mol, data['seq'], data['chirality'])
        total_chiral += sv['n_chiral']
        total_viol += sv['violations']
        nd = sum(1 for c in data['chirality'] if c == 'D')
        nl = sum(1 for c in data['chirality'] if c == 'L')
        print(f"  {sid:<16} D={nd:>2} L={nl:>2} viol={sv['violations']}")

    rate = total_viol / max(total_chiral, 1)
    print(f"\nTotal: {total_chiral} chiral residues, "
          f"{total_viol} violations ({rate:.2%})")


if __name__ == '__main__':
    main()
