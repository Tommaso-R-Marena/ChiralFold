# Changelog

## Unreleased

### Added
- **RCSB + AlphaFold DB fetch** (`chiralfold.fetch`): PDB ID, UniProt accession, or `AF-…-Fn` → local PDB for audit/correct/mirror (CLI, Web, HF Space when online).
- **mmCIF + FASTA inputs**: core `structure_io` converts mmCIF→PDB (no gemmi); FASTA with a UniProt header resolves via AFDB.
- CLI: `chiralfold fetch`, `audit --id`, `correct-af3 --id`; structure args accept `.pdb/.cif/.fasta`.

### Changed
- Overleaf/Bioinformatics package synced: clashscore + RCSB/AFDB fetch + Light/Dark Space notes in main Methods/Availability; SI Methods for clashscore and structure I/O; cover letter dated; highlights updated; zip rebuilt (main ~7 pp, SI ~10 pp).

### Fixed
- **Clash score false positives:** exclude covalent 1-2/1-3/1-4 via amino-acid topology (not a brittle 2.6 Å cutoff), skip proline amide H, ignore disulfides and donor–acceptor H-bonds, and fix amide-H placement. AFDB/PDB audits no longer report hundreds of fake clashes (e.g. LEU CA–CG).
- Clashscore: strip deposited hydrogens (re-add backbone HN only) and read **first MODEL only** — NMR/high-res files with explicit H no longer score every C–H bond as a clash.
- Regenerated `results/molprobity_comparison.json`, paper data copy, and `results/af3_experimental_systems.json` after the clashscore fix (mean `cf_clash` ~265 → ~24; panel wwPDB mean ~18).
- Windows CI `UnicodeDecodeError` when reading Colab notebooks (`tests/test_colab_publication.py` now forces UTF-8).
- `conda-recipe/meta.yaml` version bumped **3.4.0 → 3.5.1** to match PyPI.

### Changed
- README covers all headline results, Lean 4 proofs, and setup avenues (pip, conda, Docker, HF Space, Colab).
- Removed tracked bulk PDB caches (`results/d_survey/*.pdb`, `results/pdb50/*.pdb`), draft `ChiralFold_3/`, and duplicate root assets.
- Full live **mmCIF-only D-residue universe survey** (`benchmarks/mmcif_d_residue_expansion.py --mode both`): discovers RCSB entries lacking legacy `.pdb`, scans with gemmi; recovers 29/29 known errors and reports **9BC4** (DLE) as a new clear mismatch.
- `results/mmcif_only_universe_ids.json` discovery artefact.


All notable changes to this project are documented in this file.

## [Unreleased]

### Added
- mmCIF Colab: `demos/Reproduce_mmCIF_D_Residue_Survey.ipynb` (16 known-error structures → 29/29)
- Publication Colab schema CI (`tests/test_colab_publication.py`)
- Completed Aristotle Lean 4 package: orthogonal generalization, Examples.lean,
  paper Methods snippet (`formal/chirality_nogo/`)
- `docs/PYPI_PUBLISHING.md` — Trusted Publisher claim checklist + API-token fallback

### Fixed
- Results Dashboard AF3 cell (`correction.violations_after`); Lean “in progress” text
- Colab installs prefer `pip install chiralfold==3.5.1`
- Bioinformatics submission package bumped to v3.5.1; mmCIF limitation closed
- Gradio 4 / Python 3.9 CI: remove unsupported `gr.HTML(padding=...)` kwarg
- PyPI publish workflows accept optional `PYPI_API_TOKEN` and fail clearly when
  Trusted Publisher is not linked (`invalid-publisher`)

## [3.5.1] - 2026-07-20

### Added
- Reviewer-ready install path, 5-minute D-residue reproduce, clash-safety demos
- Hugging Face Space web UI with high-contrast theme + E2E web tests
- Cross-platform CI (Linux / macOS / Windows)

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
