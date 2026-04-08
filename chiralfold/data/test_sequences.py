"""
ChiralFold Test Sequence Library
=================================

30 pure D-peptide sequences + 15 mixed L/D diastereomer sequences.

Diastereomer test cases are designed to mirror the Childs, Zhou & Donald (2025)
partial-chirality benchmark scenarios where AF3 showed 50–52% violation rates.
"""

# ═══════════════════════════════════════════════════════════════════════════
# 30 Pure D-Peptide Sequences (original ChiralFold v1 suite)
# ═══════════════════════════════════════════════════════════════════════════

PURE_D_SEQS = {
    # Short (3–5 residues)
    'S01': 'FWK',
    'S02': 'RVED',
    'S03': 'ACMF',
    'S04': 'NQHWY',
    'S05': 'LFMAE',
    # Medium (6–10 residues)
    'M01': 'AFWKLD',
    'M02': 'RVFDSN',
    'M03': 'MWKYELR',
    'M04': 'FHCDSWTK',
    'M05': 'ETFSDLWKLL',
    'M06': 'TNWYQGLRFD',
    'M07': 'PRFWEYNLKA',
    'M08': 'DLKWFATINR',
    # Long (12–18 residues)
    'L01': 'FYWKELDRSNTQ',
    'L02': 'AWVELDKFRSHTN',
    'L03': 'RFSDELWNKYAMQ',
    'L04': 'THWKFVELRDSNYQA',
    'L05': 'TSFAEYWNLLSPRKD',
    'L06': 'DWFKELAYNSRTMHQV',
    'L07': 'PRVFWEYNLKASDTQMC',
    'L08': 'THWKFVELRDSNYQAMCI',
    # Homopolymers
    'H01': 'AAAAAAAAAA',
    'H02': 'FFFFFFFF',
    'H03': 'LLLLLLLLLL',
    # Charged
    'C01': 'KKRRDDEEHHNN',
    'C02': 'SSTTCCMMWWYY',
    # All-20
    'A01': 'ACDEFHIKLMNQRSTVWY',
    # Proline-rich
    'P01': 'PPVFWAEL',
    'P02': 'RFPSDPLWKP',
    # Glycine-rich (fewer chiral centers)
    'G01': 'GGGFWKGGG',
    'G02': 'AGFGWGKGL',
}

# ═══════════════════════════════════════════════════════════════════════════
# 15 Mixed L/D Diastereomer Sequences (ChiralFold v2 extension)
# ═══════════════════════════════════════════════════════════════════════════
#
# These test cases are inspired by the Childs et al. (2025) study:
#   - DP19:L-19437  (19 residues, AF3 avg 52% violation)
#   - DP9:Streptavidin (9 residues, AF3 avg 51% violation)
#   - DP12:MDM2 (12 residues, AF3 avg 50% violation)
#
# AF3 treats every residue as ~coin-flip for chirality, regardless of the
# specified L/D pattern.  ChiralFold encodes chirality explicitly.

DIASTEREOMER_SEQS = {
    # ── Childs 2025 system analogues ──────────────────────────────────────
    # DP19-inspired: 19-mer with mixed chirality
    'DIA_DP19': {
        'seq':       'ETFSDLWKLLRFNQHWYAM',
        'chirality': 'DDLLDDDLLDDDLLDDDDL',
        'note': 'DP19:L-19437 analogue (AF3: 52% violation)',
    },
    # DP9-inspired: 9-mer with alternating chirality
    'DIA_DP9': {
        'seq':       'TNWYQGLRF',
        'chirality': 'DLDLDLDLD',
        'note': 'DP9:Streptavidin analogue (AF3: 51% violation)',
    },
    # DP12-inspired: 12-mer with block chirality
    'DIA_DP12': {
        'seq':       'ETFSDLWKLLRF',
        'chirality': 'DDDDLLLLLDDL',
        'note': 'DP12:MDM2 analogue (AF3: 50% violation)',
    },

    # ── Alternating L/D patterns ──────────────────────────────────────────
    'ALT_6': {
        'seq':       'AFWKEL',
        'chirality': 'DLDLDL',
        'note': 'Alternating 6-mer',
    },
    'ALT_10': {
        'seq':       'ETFSDLWKLL',
        'chirality': 'DLDLDLDLDL',
        'note': 'Alternating 10-mer',
    },
    'ALT_14': {
        'seq':       'THWKFVELRDSNYQ',
        'chirality': 'DLDLDLDLDLDLDL',
        'note': 'Alternating 14-mer',
    },

    # ── Random L/D patterns ───────────────────────────────────────────────
    'RND_8': {
        'seq':       'MWKYELRF',
        'chirality': 'DDDLLDLD',
        'note': 'Random 8-mer',
    },
    'RND_12': {
        'seq':       'RFSDELWNKYAM',
        'chirality': 'DLDDLDDLDLDD',
        'note': 'Random 12-mer',
    },
    'RND_16': {
        'seq':       'DWFKELAYNSRTMHQV',
        'chirality': 'DLDDDLLDDLDDLDLD',
        'note': 'Random 16-mer',
    },

    # ── Drug-design scenarios (mostly L, a few D substitutions) ───────────
    'DRUG_L1D': {
        'seq':       'AFWKELDR',
        'chirality': 'LLLDLLLL',
        'note': 'Single D-substitution (drug design)',
    },
    'DRUG_L2D': {
        'seq':       'TNWYQGLRFD',
        'chirality': 'LLLDLLLLDL',
        'note': 'Two D-substitutions',
    },
    'DRUG_L3D': {
        'seq':       'RFSDELWNKYAMQ',
        'chirality': 'LLDLLLLLLDLLL',
        'note': 'Three D-substitutions',
    },

    # ── Block chirality patterns ──────────────────────────────────────────
    'BLK_HALF': {
        'seq':       'ETFSDLWKLL',
        'chirality': 'DDDDDLLLLL',
        'note': 'Half-D / half-L block',
    },
    'BLK_THIRD': {
        'seq':       'FYWKELDRSNTQ',
        'chirality': 'DDDDLLLLDDDD',
        'note': 'Thirds pattern',
    },

    # ── All-L control ─────────────────────────────────────────────────────
    'CTRL_ALL_L': {
        'seq':       'AFWKELDR',
        'chirality': 'LLLLLLLL',
        'note': 'All-L control (native protein chirality)',
    },
}

# ═══════════════════════════════════════════════════════════════════════════
# AF3 reference violation rates from Childs et al. (2025)
# ═══════════════════════════════════════════════════════════════════════════

AF3_REFERENCE = {
    'DP19:L-19437': {
        'n_residues': 19,
        'avg_violation': 0.52,
        'best_sample': 0.37,
        'worst_sample': 0.63,
        'n_experiments': 1280,
    },
    'DP9:Streptavidin': {
        'n_residues': 9,
        'avg_violation': 0.51,
        'best_sample': 0.33,
        'worst_sample': 0.67,
        'n_experiments': 1280,
    },
    'DP12:MDM2': {
        'n_residues': 12,
        'avg_violation': 0.50,
        'best_sample': 0.17,
        'worst_sample': 0.83,
        'n_experiments': 695,
    },
    'Apo_D-SH3': {
        'n_residues': 57,
        'avg_violation': 0.44,
        'best_sample': None,
        'worst_sample': None,
        'n_experiments': None,
    },
    'Synthetic_Ub_GB1': {
        'n_residues': None,
        'avg_violation': 0.44,
        'best_sample': None,
        'worst_sample': None,
        'n_experiments': None,
    },
    'Fluorinated': {
        'n_residues': None,
        'avg_violation': 0.33,
        'best_sample': None,
        'worst_sample': None,
        'n_experiments': None,
    },
    'Population_average': {
        'n_residues': None,
        'avg_violation': 0.503,
        'std_violation': 0.267,
        'n_experiments': 3195,
        'note': 'No error-free ligand in any of 100% of samples',
    },
}
