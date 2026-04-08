"""
ChiralFold — Diastereomer Enumeration
=======================================

Enumerates L/D assignments across stereogenic positions in a peptide,
generates conformers for each, audits quality, and ranks by score.

Scoring heuristic
-----------------
Diastereomers are scored with a bias toward mostly-L patterns (drug-like)
while preserving diversity of D-containing candidates:

  score = 100 - (n_d / n_stereo * 10) + conformer_count_bonus

This means pure-L scores ~100, and adding D residues mildly penalises
while more generated conformers reward. Useful for drug design screens.
"""

from __future__ import annotations

import itertools
import random
import warnings
from typing import Dict, List, Optional, Tuple

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# Maximum exhaustive enumeration limit (2^n_stereo must be ≤ this)
MAX_EXHAUSTIVE = 1000

# For long peptides (>15 residues), limit random sampling
RANDOM_SAMPLE_LIMIT = 200

# Conformer count bonus per confirmed generated conformer
CONFORMER_BONUS_PER = 0.5


# ─────────────────────────────────────────────────────────────────────────────
# Pattern generation
# ─────────────────────────────────────────────────────────────────────────────

def _stereo_positions(sequence: str) -> List[int]:
    """
    Return indices of stereogenic positions (all non-Glycine residues).

    Glycine is achiral (no Cβ), so it is excluded from the enumeration.
    Both L and D forms of all other standard amino acids are valid.
    """
    return [i for i, aa in enumerate(sequence) if aa.upper() != "G"]


def _all_patterns(n_stereo: int) -> List[Tuple[str, ...]]:
    """Enumerate all 2^n_stereo L/D patterns exhaustively."""
    return list(itertools.product("LD", repeat=n_stereo))


def _random_patterns(n_stereo: int, n_samples: int,
                     rng: random.Random) -> List[Tuple[str, ...]]:
    """
    Sample *n_samples* unique L/D patterns uniformly at random.

    If 2^n_stereo < n_samples, all patterns are returned (exhaustive).
    """
    total = 2 ** n_stereo
    if total <= n_samples:
        return _all_patterns(n_stereo)

    seen = set()
    patterns = []
    attempts = 0
    max_attempts = n_samples * 20

    while len(patterns) < n_samples and attempts < max_attempts:
        pat = tuple(rng.choice("LD") for _ in range(n_stereo))
        if pat not in seen:
            seen.add(pat)
            patterns.append(pat)
        attempts += 1

    return patterns


def _expand_pattern(
    sequence: str,
    stereo_positions: List[int],
    pattern: Tuple[str, ...],
) -> str:
    """
    Expand a compact pattern over stereo_positions into a full per-residue
    chirality string of length len(sequence).

    Glycine positions get 'L' (treated as L even though achiral — irrelevant).
    """
    full = ["L"] * len(sequence)
    for pos, chir in zip(stereo_positions, pattern):
        full[pos] = chir
    return "".join(full)


# ─────────────────────────────────────────────────────────────────────────────
# Scoring
# ─────────────────────────────────────────────────────────────────────────────

def _score_diastereomer(
    n_d: int,
    n_stereo: int,
    n_conformers_generated: int,
) -> float:
    """
    Compute a diastereomer quality score.

    Higher is better.

    Formula:
      score = 100 - (n_d / n_stereo * 10) + conformer_bonus

    where conformer_bonus = n_conformers_generated * CONFORMER_BONUS_PER

    Properties:
      - Pure-L (n_d=0): starts at 100 + conformer_bonus
      - Pure-D (n_d=n_stereo): starts at 90 + conformer_bonus
      - More generated conformers (broader local minima) → higher score
    """
    if n_stereo == 0:
        return 100.0 + n_conformers_generated * CONFORMER_BONUS_PER
    d_penalty     = (n_d / n_stereo) * 10.0
    conformer_bonus = n_conformers_generated * CONFORMER_BONUS_PER
    return round(100.0 - d_penalty + conformer_bonus, 3)


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def enumerate_diastereomers(
    sequence: str,
    top_n: int = 10,
    n_conformers: int = 3,
    seed: int = 42,
) -> List[Dict]:
    """
    Enumerate L/D diastereomers of *sequence*, generate conformers, and rank.

    Strategy:
      - Find all non-Glycine positions (stereogenic).
      - If 2^n_stereo ≤ 1000: exhaustive enumeration.
      - If 2^n_stereo > 1000 or sequence length > 15: random sampling of 200.
      - For each pattern: build SMILES, parse with RDKit, generate conformers.
      - Score using _score_diastereomer heuristic.
      - Return top_n results sorted by score descending.

    Args:
        sequence:    One-letter amino acid sequence (e.g. 'AFWKELDR').
        top_n:       Number of top-ranked diastereomers to return.
        n_conformers: Number of conformers to generate per diastereomer.
        seed:        Random seed for reproducible sampling.

    Returns:
        List of dicts (length ≤ top_n), sorted by score descending. Each dict:
          - rank (int): 1-indexed rank
          - chirality_pattern (str): Full per-residue 'L'/'D' string
          - n_d (int): Number of D residues
          - n_l (int): Number of L residues
          - n_stereo (int): Number of stereogenic positions
          - score (float): Composite score
          - smiles (str | None): SMILES string (None if RDKit parse failed)
          - n_conformers_generated (int): Conformers successfully generated
          - valid (bool): True if SMILES parsed and ≥1 conformer generated
    """
    from chiralfold.model import mixed_peptide_smiles, ChiralFold

    sequence = sequence.upper()
    stereo_positions = _stereo_positions(sequence)
    n_stereo = len(stereo_positions)

    rng = random.Random(seed)

    # Decide exhaustive vs. random sampling
    use_random = (2 ** n_stereo > MAX_EXHAUSTIVE) or (len(sequence) > 15)

    if use_random:
        sample_limit = RANDOM_SAMPLE_LIMIT
        compact_patterns = _random_patterns(n_stereo, sample_limit, rng)
    else:
        compact_patterns = _all_patterns(n_stereo)

    model = ChiralFold(n_conformers=n_conformers, fix_planarity=True)

    results = []

    for compact_pat in compact_patterns:
        full_pattern = _expand_pattern(sequence, stereo_positions, compact_pat)

        n_d = sum(1 for c in full_pattern if c == "D")
        n_l = sum(1 for c in full_pattern if c == "L")

        # Build SMILES
        smiles = None
        valid_smiles = False
        n_conf_generated = 0

        try:
            smiles = mixed_peptide_smiles(sequence, full_pattern)
            valid_smiles = True
        except Exception:
            pass

        # Generate conformers
        if valid_smiles and smiles:
            try:
                pred_result = model.predict(sequence, chirality_pattern=full_pattern)
                n_conf_generated = pred_result.get("n_conformers", 0)
            except Exception:
                n_conf_generated = 0

        score = _score_diastereomer(n_d, n_stereo, n_conf_generated)

        results.append({
            "chirality_pattern":    full_pattern,
            "n_d":                  n_d,
            "n_l":                  n_l,
            "n_stereo":             n_stereo,
            "score":                score,
            "smiles":               smiles,
            "n_conformers_generated": n_conf_generated,
            "valid":                valid_smiles,
        })

    # Sort by score descending, then fewest D residues as tiebreaker
    results.sort(key=lambda d: (-d["score"], d["n_d"]))

    # Assign ranks and return top_n
    top = results[:top_n]
    for i, r in enumerate(top, 1):
        r["rank"] = i

    return top


def format_enumeration_results(results: List[Dict]) -> str:
    """
    Pretty-print the diastereomer enumeration ranking table.

    Args:
        results: Output from enumerate_diastereomers().

    Returns:
        Multi-line formatted table string.
    """
    if not results:
        return "(no results)"

    header = (
        f"{'Rank':<5} {'Pattern':<35} {'D':>4} {'L':>4} "
        f"{'Score':>8} {'Confs':>6} {'Valid':<6}"
    )
    sep = "-" * len(header)
    lines = [
        "ChiralFold — Diastereomer Enumeration Results",
        "=" * len(header),
        header,
        sep,
    ]

    for r in results:
        pattern_str = r["chirality_pattern"]
        if len(pattern_str) > 33:
            pattern_str = pattern_str[:30] + "..."

        valid_str = "yes" if r.get("valid") else "no"

        lines.append(
            f"{r.get('rank', '?'):<5} "
            f"{pattern_str:<35} "
            f"{r['n_d']:>4d} "
            f"{r['n_l']:>4d} "
            f"{r['score']:>8.2f} "
            f"{r.get('n_conformers_generated', 0):>6d} "
            f"{valid_str:<6}"
        )

    lines.append(sep)
    lines.append(f"Showing {len(results)} diastereomers")

    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# CLI / test block
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    import os

    # Ensure the package root is on sys.path so 'chiralfold' is importable.
    # The chiralfold/threading.py shadow of stdlib threading is a known repo
    # issue; running as a module (python -m chiralfold.enumerate) avoids it.
    _pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _pkg_root not in sys.path:
        sys.path.insert(0, _pkg_root)

    # Default test: AFWK (4 residues, 4 stereogenic = 16 patterns, exhaustive)
    seq = sys.argv[1] if len(sys.argv) > 1 else "AFWK"
    top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    print("=" * 70)
    print(f"ChiralFold — Diastereomer Enumeration Self-Test")
    print(f"Sequence: {seq}  (length {len(seq)})")
    print("=" * 70)

    stereo_pos = _stereo_positions(seq)
    n_stereo   = len(stereo_pos)
    n_total    = 2 ** n_stereo
    use_random = n_total > MAX_EXHAUSTIVE or len(seq) > 15

    print(f"\nStereogenic positions: {stereo_pos} ({n_stereo} total)")
    print(f"Total possible patterns: {n_total}")
    print(f"Strategy: {'random sampling (200)' if use_random else 'exhaustive'}")
    print(f"Generating top {top_n} diastereomers ...\n")

    results = enumerate_diastereomers(seq, top_n=top_n, n_conformers=3)

    print(format_enumeration_results(results))

    if results:
        best = results[0]
        print(f"\nBest diastereomer:")
        print(f"  Pattern:   {best['chirality_pattern']}")
        print(f"  D residues: {best['n_d']}")
        print(f"  L residues: {best['n_l']}")
        print(f"  Score:      {best['score']}")
        if best.get("smiles"):
            smiles = best["smiles"]
            display = smiles if len(smiles) <= 80 else smiles[:77] + "..."
            print(f"  SMILES:    {display}")
