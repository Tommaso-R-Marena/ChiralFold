# Changelog

All notable changes to this project are documented in this file.

## [3.4.0] - 2026-06-22

### Added
- Docker image (`Dockerfile`, `docker-compose.yml`) for reproducible CLI usage
- Lean 4 formal chirality no-go proofs (`formal/chirality_nogo/`)
- Era-representative Ramachandran benchmark (n=362, chain-ID-fixed mmCIF converter)
- Full journal manuscript and bioRxiv preprint with expanded statistics
- GitHub Actions release workflow (PyPI sdist/wheel + GHCR Docker image)
- Pre-commit configuration, hardened CI coverage gate, and conda environment/recipe scaffolding
- One-click Bioinformatics Overleaf package plus fast toy dataset Colab demo
- AF3-mimetic resource benchmark documenting 100% correction recall, 0% residual violations, and ~37 ms total correction pipeline time

### Fixed
- mmCIF→PDB chain-ID truncation bug for multi-character chain identifiers
- Non-protein outlier exclusion documented for benchmark (5M2K)
- Package metadata and lock-file notes aligned to v3.4.0

### Changed
- Ramachandran headline benchmark: Spearman ρ=0.52, Pearson r=0.853 (n=362)
- Paper and preprint updated with formal verification and expanded MolProbity calibration
- README reorganized for install, API/CLI, demos, PyMOL visualization, testing, and submission-readiness workflows

## [3.2.3] - 2026-04-17

Prior stable release with benchmark reproducibility artifacts.
