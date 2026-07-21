# Docs index

| Doc | Purpose |
|-----|---------|
| [INSTALL.md](INSTALL.md) | Cross-platform pip / conda / Docker install |
| [../README.md](../README.md) | Project overview + reviewer path |
| [../results/DATA_README.md](../results/DATA_README.md) | Results artefacts |
| [../results/REPRODUCIBILITY.md](../results/REPRODUCIBILITY.md) | Benchmark reproduction |
| [../benchmarks/README.md](../benchmarks/README.md) | Benchmark scripts |
| [../paper/submission/bioinformatics/](../paper/submission/bioinformatics/) | Manuscript package |
| [REVIEWER_COLAB.md](REVIEWER_COLAB.md) | **Reviewer path:** Colab notebooks + HF Space (ordered) |
| [BRANCH_PROTECTION.md](BRANCH_PROTECTION.md) | Require **CI success** before merging to master |

## Colab notebooks (`demos/`)

| Notebook | Time | Purpose |
|----------|------|---------|
| `Reproduce_PDB_D_Residue_Errors_5min.ipynb` | &lt;5 min | Offline survey recompute (12,573 / 29 errors) |
| `Demo_Unusual_Cases_Clash_Safety.ipynb` | ~3 min | Macrocycles, CCD ligands, clash isometry |
| `ChiralFold_Quick_Demo.ipynb` | ~2 min | Audit / predict / mirror |
| `ChiralFold_Toy_Dataset_Demo.ipynb` | ~2 min | Offline toy + AF3 fix |
| `ChiralFold_Results_Dashboard.ipynb` | ~2 min | Interactive Plotly dashboard |
| `D_Residue_Experimental_Validation.ipynb` | ~2 min | RCSB/CCD cross-check (network) |
| `Expanded_Ramachandran_Benchmark.ipynb` | longer | Full Rama expansion |

> The two reviewer-facing notebooks have **distinct filenames** so Colab downloads are not confused (`Reproduce_…` vs `Demo_Unusual_…`).
