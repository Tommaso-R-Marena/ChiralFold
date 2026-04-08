"""
ChiralFold — Command-Line Interface
=====================================

Usage:
    chiralfold predict SEQUENCE [--chirality PATTERN]
    chiralfold validate SEQUENCE --chirality PATTERN
    chiralfold audit PDB_FILE [--chain CHAIN] [--json]
    chiralfold audit --rcsb-batch FILE -o OUTPUT.csv
    chiralfold correct-af3 PDB_FILE [--output OUTPUT]
    chiralfold enumerate SEQUENCE [--top N]
    chiralfold score-interface RECEPTOR LIGAND
    chiralfold mirror PDB_FILE [--output OUTPUT] [--chains A,B]
    chiralfold mirror-id PDB_ID [--output OUTPUT]
    chiralfold benchmark
"""

import argparse
import csv
import json
import os
import sys

from . import __version__
from .model import ChiralFold, mixed_peptide_smiles
from .validator import validate_diastereomer


def main():
    parser = argparse.ArgumentParser(
        prog='chiralfold',
        description=(
            'ChiralFold: Stereochemistry toolkit for proteins and peptides'
        ),
    )
    parser.add_argument(
        '--version', action='version', version=f'chiralfold {__version__}'
    )

    sub = parser.add_subparsers(dest='command', help='Available commands')

    # ── predict ───────────────────────────────────────────────────────────
    p = sub.add_parser('predict', help='Predict peptide structure')
    p.add_argument('sequence', help='Amino acid sequence')
    p.add_argument('--chirality', '-c', default=None)
    p.add_argument('--json', action='store_true')

    # ── validate ──────────────────────────────────────────────────────────
    p = sub.add_parser('validate', help='Validate peptide chirality')
    p.add_argument('sequence')
    p.add_argument('--chirality', '-c', required=True)
    p.add_argument('--json', action='store_true')

    # ── audit ─────────────────────────────────────────────────────────────
    p = sub.add_parser('audit', help='Audit PDB structure quality')
    p.add_argument('pdb_file', nargs='?', default=None, help='PDB file path')
    p.add_argument('--chain', default=None)
    p.add_argument('--json', action='store_true')
    p.add_argument('--rcsb-batch', dest='rcsb_batch', default=None,
                   help='File with PDB IDs (one per line) for batch RCSB audit')
    p.add_argument('-o', '--output', default=None,
                   help='Output CSV for batch mode')

    # ── correct-af3 ──────────────────────────────────────────────────────
    p = sub.add_parser('correct-af3',
                       help='Correct chirality violations in AF3 output')
    p.add_argument('pdb_file', help='AF3 prediction PDB file')
    p.add_argument('--output', '-o', default=None)

    # ── enumerate ─────────────────────────────────────────────────────────
    p = sub.add_parser('enumerate',
                       help='Enumerate and rank diastereomers')
    p.add_argument('sequence', help='Amino acid sequence')
    p.add_argument('--top', type=int, default=10)
    p.add_argument('--json', action='store_true')

    # ── score-interface ───────────────────────────────────────────────────
    p = sub.add_parser('score-interface',
                       help='Score binding interface')
    p.add_argument('receptor', help='Receptor PDB file')
    p.add_argument('ligand', help='Ligand PDB file')
    p.add_argument('--receptor-chain', default=None)
    p.add_argument('--ligand-chain', default=None)
    p.add_argument('--json', action='store_true')

    # ── mirror ────────────────────────────────────────────────────────────
    p = sub.add_parser('mirror', help='Mirror PDB structure (L↔D)')
    p.add_argument('pdb_file')
    p.add_argument('--output', '-o', default=None)
    p.add_argument('--chains', default=None)

    # ── mirror-id ─────────────────────────────────────────────────────────
    p = sub.add_parser('mirror-id', help='Download and mirror PDB by ID')
    p.add_argument('pdb_id')
    p.add_argument('--output', '-o', default=None)
    p.add_argument('--chains', default=None)

    # ── benchmark ─────────────────────────────────────────────────────────
    sub.add_parser('benchmark', help='Run full benchmark suite')

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        sys.exit(0)

    handlers = {
        'predict': _cmd_predict, 'validate': _cmd_validate,
        'audit': _cmd_audit, 'correct-af3': _cmd_correct_af3,
        'enumerate': _cmd_enumerate, 'score-interface': _cmd_score_interface,
        'mirror': _cmd_mirror, 'mirror-id': _cmd_mirror_id,
        'benchmark': _cmd_benchmark,
    }
    handlers[args.command](args)


# ═══════════════════════════════════════════════════════════════════════════
# Command handlers
# ═══════════════════════════════════════════════════════════════════════════

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


def _cmd_audit(args):
    from .auditor import audit_pdb, format_report

    # Batch mode
    if args.rcsb_batch:
        _audit_rcsb_batch(args.rcsb_batch, args.output)
        return

    if args.pdb_file is None:
        print("Error: provide a PDB file or use --rcsb-batch", file=sys.stderr)
        sys.exit(1)

    report = audit_pdb(args.pdb_file)
    if args.json:
        safe = {k: v for k, v in report.items() if k not in ('_atoms',)}
        print(json.dumps(safe, indent=2, default=str))
    else:
        print(format_report(report))


def _audit_rcsb_batch(batch_file, output_csv):
    """Batch audit PDB IDs from a file, output CSV."""
    import urllib.request
    import tempfile
    from .auditor import audit_pdb

    with open(batch_file) as f:
        pdb_ids = [line.strip().upper() for line in f
                   if line.strip() and not line.startswith('#')]

    print(f"Auditing {len(pdb_ids)} structures from RCSB...")

    rows = []
    for i, pdb_id in enumerate(pdb_ids):
        if len(pdb_id) != 4:
            continue
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
                urllib.request.urlretrieve(url, tmp.name)
                tmp_path = tmp.name

            report = audit_pdb(tmp_path)
            rows.append({
                'pdb_id': pdb_id,
                'chirality_pct': report['chirality']['pct_correct'],
                'rama_favored_pct': report['ramachandran']['pct_favored'],
                'rama_outlier_pct': report['ramachandran']['pct_outlier'],
                'planarity_pct': report['planarity']['pct_within_6deg'],
                'clash_score': report['clashes']['clash_score'],
                'overall_score': report['overall_score'],
            })
            os.unlink(tmp_path)

            if (i + 1) % 5 == 0:
                print(f"  {i+1}/{len(pdb_ids)} done...")

        except Exception as e:
            rows.append({
                'pdb_id': pdb_id, 'chirality_pct': 'ERROR',
                'rama_favored_pct': '', 'rama_outlier_pct': '',
                'planarity_pct': '', 'clash_score': '',
                'overall_score': str(e)[:40],
            })

    # Write CSV
    out_path = output_csv or 'audit_results.csv'
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'pdb_id', 'chirality_pct', 'rama_favored_pct',
            'rama_outlier_pct', 'planarity_pct', 'clash_score',
            'overall_score'
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved {len(rows)} results to {out_path}")


def _cmd_correct_af3(args):
    from .af3_correct import correct_af3_output

    output = args.output
    if output is None:
        base = args.pdb_file.rsplit('.', 1)[0]
        output = f"{base}_corrected.pdb"

    result = correct_af3_output(args.pdb_file, output)

    if 'error' in result:
        print(f"Error: {result['error']}", file=sys.stderr)
        sys.exit(1)

    n_before = result.get('violations_before', {}).get('n_violations', 0)
    n_after = result.get('violations_after', {}).get('n_violations', 0)
    print(f"Violations before: {n_before}")
    print(f"Corrected:         {result.get('n_corrected', 0)}")
    print(f"Violations after:  {n_after}")
    print(f"Output: {output}")


def _cmd_enumerate(args):
    from .enumerate import enumerate_diastereomers, format_enumeration_results

    seq = args.sequence.upper()
    results = enumerate_diastereomers(seq, top_n=args.top)

    if args.json:
        print(json.dumps(results, indent=2, default=str))
    else:
        format_enumeration_results(results)


def _cmd_score_interface(args):
    from .interface_scorer import score_interface

    result = score_interface(
        args.receptor, args.ligand,
        receptor_chain=args.receptor_chain,
        ligand_chain=args.ligand_chain,
    )

    if args.json:
        print(json.dumps(result, indent=2, default=str))
    else:
        bsa = result.get('bsa', result.get('buried_surface_area', 0))
        sc = result.get('shape_complementarity', {})
        sc_val = sc.get('sc_fraction', 0) if isinstance(sc, dict) else sc
        hb = result.get('hbonds', result.get('n_hbonds', 0))
        sb = result.get('salt_bridges', result.get('n_salt_bridges', 0))
        hc = result.get('hydrophobic_contacts', result.get('n_hydrophobic', 0))
        n_pairs = result.get('n_interface_pairs', 0)
        score = result.get('interface_score', 0)
        print(f"Buried Surface Area:    {bsa:.0f} A^2")
        print(f"Shape Complementarity:  {sc_val:.3f}")
        print(f"Hydrogen bonds:         {hb}")
        print(f"Salt bridges:           {sb}")
        print(f"Hydrophobic contacts:   {hc}")
        print(f"Interface pairs:        {n_pairs}")
        print(f"Interface score:        {score:.1f}/100")


def _cmd_mirror(args):
    from .pdb_pipeline import mirror_pdb
    output = args.output
    if output is None:
        base = args.pdb_file.rsplit('.', 1)[0]
        output = f"{base}_mirror.pdb"
    chains = args.chains.split(',') if args.chains else None
    result = mirror_pdb(args.pdb_file, output, chains=chains)
    print(f"Mirrored {result['n_atoms']} atoms, "
          f"{result['n_residues']} residues → {result['output_path']}")


def _cmd_mirror_id(args):
    from .pdb_pipeline import fetch_and_mirror
    output = args.output or f"{args.pdb_id.upper()}_mirror.pdb"
    chains = args.chains.split(',') if args.chains else None
    result = fetch_and_mirror(args.pdb_id, output, chains=chains)
    print(f"Downloaded {args.pdb_id.upper()}, mirrored {result['n_atoms']} atoms "
          f"→ {result['output_path']}")


def _cmd_benchmark(args):
    from .data.test_sequences import PURE_D_SEQS, DIASTEREOMER_SEQS
    from .validator import validate_smiles_chirality
    from .model import d_peptide_smiles
    from rdkit import Chem

    total_chiral = total_viol = 0
    print("Pure D-peptides:")
    for sid, seq in PURE_D_SEQS.items():
        smi = d_peptide_smiles(seq)
        mol = Chem.MolFromSmiles(smi)
        sv = validate_smiles_chirality(mol, seq, 'D' * len(seq))
        total_chiral += sv['n_chiral']
        total_viol += sv['violations']

    print("Diastereomers:")
    for sid, data in DIASTEREOMER_SEQS.items():
        smi = mixed_peptide_smiles(data['seq'], data['chirality'])
        mol = Chem.MolFromSmiles(smi)
        sv = validate_smiles_chirality(mol, data['seq'], data['chirality'])
        total_chiral += sv['n_chiral']
        total_viol += sv['violations']

    rate = total_viol / max(total_chiral, 1)
    print(f"\nTotal: {total_chiral} chiral residues, "
          f"{total_viol} violations ({rate:.2%})")


if __name__ == '__main__':
    main()
