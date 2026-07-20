# Overleaf one-click compile

## Upload

1. Download `ChiralFold_Bioinformatics_Overleaf.zip` from this folder (or from `Overleaf/` at repo root).
2. Overleaf → **New Project** → **Upload Project** → select the zip.
3. Set **Main document** to `chiralfold_bioinformatics.tex` (Menu → Main document).
4. Compiler: **pdfLaTeX** (default). Bibliography: **BibTeX**.
5. Click **Recompile**.

Supplementary: compile `supplementary_material.tex` as a second document (or upload separately).

## Contents

- `chiralfold_bioinformatics.tex` — main paper (≤7 pages target)
- `supplementary_material.tex` — methods + tables + extended figures
- `figures/fig1`–`fig7` PNGs
- `references.bib`, OUP class/bst
- Bundled `data/` for figure regeneration

## Local rebuild of the zip

```bash
cd paper/submission/bioinformatics
./package_overleaf.sh
cp ChiralFold_Bioinformatics_Overleaf.zip ../../../Overleaf/
```
