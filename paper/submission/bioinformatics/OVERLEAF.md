# Overleaf one-click compile

## Upload

1. Download `ChiralFold_Bioinformatics_Overleaf.zip` from this folder (or from `Overleaf/` at repo root).
2. Overleaf → **New Project** → **Upload Project** → select the zip.
3. Set **Main document** to `chiralfold_bioinformatics.tex` (Menu → Main document).
4. Compiler: **pdfLaTeX** (default). Bibliography: **BibTeX**.
5. Click **Recompile**.

Supplementary: compile `supplementary_material.tex` as a second document (or upload separately).

## Known harmless warnings

- **Font shape `OT1/cmss/...` undefined** — OUP template uses Computer Modern Sans sizes that TeX substitutes; cosmetic only.
- **Underfull/Overfull `\hbox` in headers** — decorative header rules in `oup-authoring-template.cls` (`modern` option); does not affect body text.
- The PNAS Nexus logo is intentionally disabled in the class file (`\societylogo{}`) so Bioinformatics compiles without that asset.

## Contents

- `chiralfold_bioinformatics.tex` — main paper (≤7 pages target; currently ~6)
- `supplementary_material.tex` — methods + Lean Notes S1 + extended figures
- `figures/fig1`–`fig7` PNGs (300 DPI)
- `references.bib`, OUP class/bst
- Bundled `data/` for figure regeneration

## Local rebuild of the zip

```bash
bash scripts/build_overleaf_package.sh
# or:
cd paper/submission/bioinformatics && ./package_overleaf.sh
```
