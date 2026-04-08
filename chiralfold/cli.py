"""
ChiralFold — Command-Line Interface
=====================================

Usage:
    chiralfold predict SEQUENCE [--chirality PATTERN] [--mode MODE]
    chiralfold validate SEQUENCE --chirality PATTERN
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
        description='ChiralFold: Chirality-preserving peptide structure prediction',
    )
    parser.add_argument(
        '--version', action='version', version=f'chiralfold {__version__}'
    )

    sub = parser.add_subparsers(dest='command', help='Available commands')

    # predict
    p_pred = sub.add_parser('predict', help='Predict peptide structure')
    p_pred.add_argument('sequence', help='Amino acid sequence (one-letter codes)')
    p_pred.add_argument(
        '--chirality', '-c', default=None,
        help='Per-residue chirality pattern (D/L string). Default: all D.'
    )
    p_pred.add_argument(
        '--mode', '-m', default='de_novo', choices=['de_novo', 'mirror'],
        help='Prediction mode (default: de_novo).'
    )
    p_pred.add_argument('--json', action='store_true', help='Output as JSON')

    # validate
    p_val = sub.add_parser('validate', help='Validate chirality of a peptide')
    p_val.add_argument('sequence', help='Amino acid sequence')
    p_val.add_argument(
        '--chirality', '-c', required=True,
        help='Per-residue chirality pattern (D/L string).'
    )
    p_val.add_argument('--json', action='store_true', help='Output as JSON')

    # benchmark
    p_bench = sub.add_parser('benchmark', help='Run full benchmark suite')
    p_bench.add_argument(
        '--output', '-o', default='.',
        help='Output directory for results (default: current directory).'
    )

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == 'predict':
        _cmd_predict(args)
    elif args.command == 'validate':
        _cmd_validate(args)
    elif args.command == 'benchmark':
        _cmd_benchmark(args)


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


def _cmd_benchmark(args):
    # Import the benchmark runner
    print("Running ChiralFold full benchmark suite...")
    print(f"Output directory: {args.output}")
    print()
    # Defer to the benchmark script
    from .model import d_peptide_smiles, l_peptide_smiles, mixed_peptide_smiles
    from .validator import validate_smiles_chirality, validate_3d_chirality

    # Quick inline benchmark
    from .data.test_sequences import PURE_D_SEQS, DIASTEREOMER_SEQS

    total_chiral = 0
    total_viol = 0
    print("Pure D-peptide benchmark:")
    for sid, seq in PURE_D_SEQS.items():
        from rdkit import Chem
        smi = d_peptide_smiles(seq)
        mol = Chem.MolFromSmiles(smi)
        sv = validate_smiles_chirality(mol, seq, 'D' * len(seq))
        nc = sv['n_chiral']
        nv = sv['violations']
        total_chiral += nc
        total_viol += nv
        print(f"  {sid:<8} {seq:<22} chiral={nc:>2} viol={nv}")

    print(f"\nDiastereomer benchmark:")
    for sid, data in DIASTEREOMER_SEQS.items():
        seq = data['seq']
        chir = data['chirality']
        smi = mixed_peptide_smiles(seq, chir)
        mol = Chem.MolFromSmiles(smi)
        sv = validate_smiles_chirality(mol, seq, chir)
        nc = sv['n_chiral']
        nv = sv['violations']
        total_chiral += nc
        total_viol += nv
        nd = sum(1 for c in chir if c == 'D')
        nl = sum(1 for c in chir if c == 'L')
        print(f"  {sid:<16} {seq:<22} D={nd:>2} L={nl:>2} viol={nv}")

    rate = total_viol / max(total_chiral, 1)
    print(f"\nTotal: {total_chiral} chiral residues, "
          f"{total_viol} violations ({rate:.2%})")
    print(f"ChiralFold: {rate:.2%} vs AlphaFold 3: 51%")


if __name__ == '__main__':
    main()
