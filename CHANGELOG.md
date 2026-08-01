# Changelog

## [3.6.0] - 2026-07-31

Correctness and performance release. No published survey number changes: the
PDB-wide D-residue result (12,573 checkable residues, 29 mismatches, 16
structures, 0.23%) is produced by `benchmarks/independent_d_residue_verification.py`,
which imports no ChiralFold code, and was re-verified after every change in this
release (`python benchmarks/reproduce_d_residue_errors.py` → PASS).

### Fixed — correctness

- **C-terminal glycine built an extra residue.** `mixed_peptide_smiles` appended
  two fragments for a trailing `G` — a malformed carbamate `NC(=O)O` *and* the
  real `NCC(=O)O` — so every sequence ending in glycine was one residue too long
  (`d_peptide_smiles("G")` returned a two-unit molecule, formula `C3H6N2O3`
  instead of glycine's `C2H5NO2`). Affects `ChiralFold.predict`,
  `enumerate_diastereomers`, and `validate_diastereomer` for any G-terminated
  sequence. The three Childs et al. benchmark sequences contain no glycine, so
  `results/childs2025_comparison.json` (478 residues, 0 violations) is unchanged.
- **chi1 sign convention was inverted.** `rotamers.compute_chi1` returned the
  *negation* of the IUPAC/BioPython dihedral, so L side chains were scored
  against effectively mirrored rotamer wells: a residue in the common g- rotamer
  was read as g+ and counted as an outlier. On 3IWY the chi1 outlier rate falls
  from 55.7% to 5.4% and favoured from 43.7% to 93.4%. `compute_chi1` now agrees
  with `auditor._dihedral_deg` to machine precision (asserted in tests).
- **D-residues had no rotamer library at all.** `validate_rotamers` reported
  `n_residues_checked = 0` for every D-peptide. D-CCD codes (DAL…DVA, MED) now
  carry the mirrored chi1 wells, and common modified residues (MSE, HIE/HID/HIP,
  CYX/CYM, SEP/TPO/PTR) alias to their parents.
- **THR chi1 library omitted its most populated well.** g+ (62°, 49% of
  threonines per Lovell 2000) was missing, so roughly half of all THR residues
  scored as chi1 outliers.
- **Cβ-less residues were classified by an invented Cβ.** `af3_correct` fell back
  to an idealised *L*-configured pseudo-Cβ whose signed volume is positive by
  construction, producing a false violation for every D-residue with an
  unmodelled side chain — and a meaningless "correction" (reflecting Cα across an
  invented plane does not invert the real stereocentre). Detection now uses Hα
  when Cβ is absent, and otherwise reports the residue as `n_unassignable`.
- **Near-planar Cα tetrahedra were forced into a class.** `af3_correct` now
  applies the same `|V| < 0.05 Å³` guard the auditor uses.
- **Violations were keyed without the insertion code.** A violation at residue 47
  also "corrected" residue 47A. Keys are now
  `(model, chain, resseq, icode)` throughout.
- **Multi-model files were mis-handled, differently, in four places.**
  `af3_correct`, `interface_scorer` and `rotamers` read every model and then
  silently collapsed them by atom-name key (an NMR ensemble was analysed as one
  structure with *n* superimposed copies); the auditor read only model 1.
  `af3_correct` now detects and corrects every model; the other readers take the
  first model explicitly.
- **Mirroring an NMR ensemble flattened it.** `mirror_pdb` dropped
  `MODEL`/`ENDMDL`, emitting an *n*-model input as a single model containing *n*
  superimposed copies of every residue, and reporting `n_atoms` inflated by *n*.
  Model framing is preserved, and `TITLE`/`CRYST1`/`COMPND` header records are
  no longer discarded.
- **Residues modelled only as altloc B/C vanished.** The per-atom accept-list
  `{' ', 'A', '1'}` deleted every atom of such residues, and where alternates did
  exist it could keep the minority conformer or mix atoms from different
  alternates within one residue. Selection is now per residue by summed
  occupancy. Scope measured in `results/altloc_policy_sensitivity.json`:
  12/13 corpus structures are bit-identical, and **Ramachandran outlier
  percentages are unchanged on every structure**, including the one with a
  majority-B conformer. `_pdbio.select_altlocs(policy="legacy")` reproduces the
  old behaviour so the effect on any number can be measured, not argued.
- **Insertion codes collided in `interface_scorer`.** De-duplication keyed on
  `(chain, resseq, name)` merged residues 47 and 47A.
- **Salt bridges were measured Cα–Cα at 4 Å**, a distance essentially never
  reached across an interface, so the term contributed zero on real complexes.
  Now measured between charged side-chain groups (Lys NZ; Arg NE/NH1/NH2;
  His ND1/NE2 vs Asp OD1/OD2, Glu OE1/OE2) per Barlow & Thornton (1983), counted
  once per residue pair, with a Cα fallback for backbone-only models.
- **Quoted mmCIF atom names kept their quotes**, so nucleic-acid names arrived as
  `"O5'` and matched nothing downstream.
- **`_random_patterns` could return fewer patterns than requested.** Rejection
  sampling gave up after `20 × n_samples` attempts; sampling is now exact.

### Changed — performance

Reproduce with `python benchmarks/compare_baseline_performance.py --baseline v3.5.1`
(spawns a git worktree at the baseline and times identical public-API code in both).

| Operation | v3.5.1 | v3.6.0 | Speedup |
|---|---|---|---|
| `audit_pdb`, 3IWY (1,589 atoms) | 87.3 ms | 22.2 ms | **3.9×** |
| `audit_pdb`, 6,356 atoms | 371.0 ms | 106.4 ms | **3.5×** |
| `audit_pdb`, 25,424 atoms | 1,628 ms | 448 ms | **3.6×** |
| `score_interface`, 1 copy | 97.6 ms | 31.6 ms | **3.1×** |
| `score_interface`, 16 copies | 1,425 ms | 233 ms | **6.1×** |
| `detect_chirality_violations` | 8.2 ms | 6.4 ms | 1.3× |
| Ramachandran classification (20k φ/ψ) | 175 ms | 6.9 ms | **25×** |
| `ChiralFold.predict` (AF3 benchmark, mean) | 17.6 s | 10.6 s | 1.7× |

- **Shared residue table.** The chirality, bond-geometry, Ramachandran and
  planarity checks rebuilt a per-residue `{name: xyz}` dict *each*; they now share
  one column-oriented `_ResidueTable` built once per audit.
- **Batched geometry kernels.** Signed volumes, bond lengths/angles and dihedrals
  are evaluated for all residues in one call each. `np.cross`/`np.linalg.norm` on
  3-vectors cost more in dispatch than arithmetic, and the audit made ~30 such
  calls per residue. Agreement with the scalar reference is asserted in tests.
- **Clash detection reordered and vectorised.** The KD-tree radius is now the
  largest distance that can actually overlap (`2·max(r) − 0.4`, not `2·max(r)`),
  cutting candidate pairs by ~27% by volume; candidates are then filtered by a
  vectorised overlap test, a packed-int `searchsorted` against the covalent
  exclusion set, and only then the per-pair H-bond geometry test (safe to defer —
  the H-bond rule can only *suppress* a clash).
- **Memoised residue topology.** Same-residue 1-2/1-3/1-4 exclusions are computed
  once per distinct `(resname, atom-name signature)` instead of by breadth-first
  search per residue.
- **Batched hydrogen placement.** Reduce-style HN/HA/HA2/HA3 placement is planned
  in one pass and evaluated in three batched calls.
- **KD-trees in `interface_scorer`.** The dense `(N_rec, N_lig)` distance matrix
  is gone: two 8,000-atom chains needed 1,349 MB and 5.6 s, and larger pairs could
  not be scored at all; the same query now takes 192 MB and 0.9 s.
- **Vectorised Ramachandran classification** (`ramachandran.classify_regions`),
  including the empirical-grid lookup.
- **Threaded force-field minimisation.** `MMFFOptimizeMoleculeConfs` /
  `UFFOptimizeMoleculeConfs` now run with `numThreads=0`; embedding was already
  multithreaded but minimisation dominated ensemble wall time.
- **Fewer parses per correction.** `correct_af3_output` parsed its input up to
  four times; it now parses once and re-reads only the file it wrote, to verify.

### Added

- `chiralfold/_pdbio.py` — one canonical fixed-column PDB reader shared by the
  auditor, corrector, interface scorer and rotamer validator, with explicit
  policies for models, water, insertion codes and alternate locations.
- `benchmarks/performance_benchmark.py` — throughput/scaling harness (offline;
  tiles a real deposited structure) → `results/performance_benchmark.json`.
- `benchmarks/compare_baseline_performance.py` — before/after comparison against
  any git revision via a throwaway worktree → `results/performance_comparison.json`.
- `benchmarks/altloc_policy_sensitivity.py` — measures the audit-level effect of
  the alternate-location policy change → `results/altloc_policy_sensitivity.json`.
- `tests/test_regressions_v360.py` — 38 regression tests, one per defect above,
  plus scalar-vs-batched equivalence for every vectorised kernel.
- `auditor.backbone_dihedrals`, `auditor.CHIRALITY_PLANAR_EPS`,
  `auditor.CLASH_OVERLAP_CUTOFF`, `ramachandran.classify_regions`, and
  `model._max_pairwise_distance_error` as documented API.
- `predict_from_mirror` / `verify_mirror_chirality` now report
  `max_pairwise_distance_error` — the *checkable* content of the isometry claim,
  where the pre-existing `rmsd_to_ideal_mirror` compares a reflection with itself
  and is 0 by construction.

### Known limitations

- The wwPDB Ramachandran calibration (*n* = 362) and the 48-structure MolProbity
  head-to-head download their inputs from RCSB and were not re-run in this
  release (no network in the release environment). Their inputs are not cached in
  the repository. The φ/ψ computation is bit-identical to v3.5.1, and the
  alternate-location change is shown above to leave Ramachandran outlier
  percentages unchanged across the repository corpus, so no change is expected;
  re-run `benchmarks/expand_ramachandran_benchmark.py` and
  `benchmarks/molprobity_comparison.py` with network access to confirm.

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
- Amide-H builder is chain- and continuity-aware (no cross-chain H from another chain's carbonyl; insertion codes in residue keys); peptide 1–4 list includes C(i)–CB(i+1).
- **Probe-style Cα hydrogens:** place HA (residues with CB) and Gly HA2/HA3; exclude peptide/same-residue 1–5 contacts involving model H; H VDW 1.17 Å.
- **Geometry-gated H-bonds:** polar H···Acc only if ≤2.5 Å and angle ≥120°; heavy N···O only in 2.5–3.5 Å; HA never treated as donor.
- Regenerated `results/molprobity_comparison.json`, paper data copy, and `results/af3_experimental_systems.json` after the clashscore fix (mean `cf_clash` ~265 → ~27; panel wwPDB mean ~18).
- Windows CI `UnicodeDecodeError` when reading Colab notebooks (`tests/test_colab_publication.py` now forces UTF-8).
- `conda-recipe/meta.yaml` version bumped **3.4.0 → 3.5.1** to match PyPI.

### Added
- `audit_pdb(..., chain=)` and CLI `chiralfold audit --chain`; `audit_pdb` accepts mmCIF via `ensure_pdb_path`.
- Batch RCSB audit uses `fetch_rcsb` (PDB + mmCIF fallback) instead of raw `.pdb` urllib.
- Web / HF Space audit KPIs show Clashscore alongside chirality / Rama.
- Tests for Gly HA2/HA3, HA non-donor, and H-bond angle/distance gates.

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
