# ChiralFold Toy Example Data

This directory contains tiny, package-distributed PDB examples for unit tests,
documentation snippets, and offline smoke checks.

## Files

- `toy_ubiquitin_fragment.pdb`: the first 10 residues of ubiquitin chain A from
  PDB entry 1UBQ (`MQIFVKTLTG`). It was extracted from `/tmp/1UBQ.pdb` when
  available during test setup. The file is intentionally small and contains
  only ATOM records plus minimal header remarks, so it can be used without
  network access.

The fragment is not intended as a benchmark-quality dataset; use the full
benchmark inputs under `results/` and `benchmarks/` for scientific validation.
