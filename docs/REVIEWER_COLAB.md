# Reviewer reproduction path (Colab + web)

Use this ordered path when evaluating the Bioinformatics submission. Total time: **~15 minutes** for the core offline checks; add ~5 minutes for optional network notebooks.

## 1. Five-minute survey (offline, numpy only)

**Notebook:** [`demos/Reproduce_PDB_D_Residue_Errors_5min.ipynb`](../demos/Reproduce_PDB_D_Residue_Errors_5min.ipynb)

Recomputes the headline **12,573 residues / 29 errors / 16 structures** from the frozen CSV. No ChiralFold import, no RCSB network after clone.

## 2. Interactive results dashboard (offline plots)

**Notebook:** [`demos/ChiralFold_Results_Dashboard.ipynb`](../demos/ChiralFold_Results_Dashboard.ipynb)

Single notebook covering all headline metrics:

| Panel | What it shows |
|-------|----------------|
| Ramachandran progression | Spearman ρ vs sample size (paper headline n=362) |
| Ramachandran scatter | ChiralFold vs wwPDB outlier % |
| D-residue taxonomy | 29 errors by class and CCD code |
| Experimental validation | 16 error structures + automated pass/borderline |
| AF3 correction | 100% detection / 0% residual on Childs-system peptides |
| mmCIF re-verification | 29/29 known errors recovered; universe scan (9BC4) |
| Lean / interpretation | Machine-checked no-go theorems summary |

Manifest: `results/colab_integrated_manifest.json` (version 3.5.1).

## 3. mmCIF expansion (optional; network)

**Notebook:** [`demos/Reproduce_mmCIF_D_Residue_Survey.ipynb`](../demos/Reproduce_mmCIF_D_Residue_Survey.ipynb)

Live RCSB mmCIF scan: known-error cohort + mmCIF-only universe ($n=79$ at freeze).

## 4. Experimental validation (optional; network)

**Notebook:** [`demos/D_Residue_Experimental_Validation.ipynb`](../demos/D_Residue_Experimental_Validation.ipynb)

RCSB metadata + wwPDB CCD InChI cross-check for all 16 error structures.

## 5. Edge cases + clash isometry (~3 min)

**Notebook:** [`demos/Demo_Unusual_Cases_Clash_Safety.ipynb`](../demos/Demo_Unusual_Cases_Clash_Safety.ipynb)

Macrocycles (1XT7, 2RMI), CCD ligands (1OF6), high-resolution density (1HHZ).

## 6. Web UI (no install)

**Hugging Face Space (light UI):** https://huggingface.co/spaces/The-Philosopher/ChiralFold-App?__theme=light

Use the `?__theme=light` URL so Gradio does not follow OS dark mode (that was causing unreadable file names and audit reports).

Upload PDB → Audit → Correct → Mirror. Same workflows tested in CI (`tests/test_web_e2e.py`).

## 7. Formal proofs (optional; local)

```bash
cd formal/chirality_nogo && lake build
```

Machine-checked Lean 4 no-go theorems; full derivations in Supplementary Notes S1.

## Scope reminder for reviewers

ChiralFold is a **stereochemistry auditor and coordinate-mapping post-processor**. It does **not** model thermodynamic ensembles, sample Boltzmann-weighted conformers, or perform molecular-mechanics energy minimization. See `formal/chirality_nogo/paper_methods_snippet.tex` (Limitation paragraph) and Discussion §Scope and limitations in the main manuscript.

## CI coverage

- `tests/test_colab_publication.py` — notebook schema + manifest freshness
- `tests/test_web_e2e.py` — web + HF Space audit/correct/mirror parity
- `.github/workflows/ci.yml` — HF Space Docker smoke on every PR
- `.github/workflows/deploy-hf-space.yml` — deploy to HF after tests (upstream repo only)
