---
title: ChiralFold
emoji: 🧬
colorFrom: green
colorTo: blue
sdk: docker
app_port: 7860
pinned: true
license: mit
tags:
  - biology
  - chemistry
  - bioinformatics
  - protein
  - stereochemistry
  - alphafold
  - pdb
short_description: Audit, correct, and mirror PDB stereochemistry
---

# ChiralFold — Stereochemistry Web App

**Live:** https://huggingface.co/spaces/The-Philosopher/ChiralFold-App

Use the in-app **Light / Dark** toggle (or `?__theme=light` / `?__theme=dark`). Both themes keep file names, audit reports, and the footer readable.

Upload any PDB structure to:

1. **Audit** Cα chirality, Ramachandran, planarity, clashes
2. **Correct** inverted stereocenters (AF3 / diffusion outputs)
3. **Mirror** for exact L↔D binder / D-peptide design

## Links

- [GitHub](https://github.com/Tommaso-R-Marena/ChiralFold)
- Landing redirect: https://huggingface.co/spaces/The-Philosopher/ChiralFold
- Local: `pip install "chiralfold[web]" && chiralfold-web`

**Privacy:** files are processed in this Space session and are not retained.

<!-- deploy-stamp: 2026-07-21-theme-toggle-v1 -->
