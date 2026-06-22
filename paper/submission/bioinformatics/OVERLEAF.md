# Overleaf upload — ChiralFold Bioinformatics submission

## Quick start

1. Run locally (optional):
   ```bash
   ./package_overleaf.sh
   ```
2. Upload `ChiralFold_Bioinformatics_Overleaf.zip` to [Overleaf](https://www.overleaf.com) → **New Project** → **Upload Project**.
3. Set **Main document** to `chiralfold_bioinformatics.tex` (Menu → Main document).
4. Click **Recompile**. Use **pdfLaTeX** + **BibTeX** (Overleaf default).

## What is included

| Path | Purpose |
|------|---------|
| `chiralfold_bioinformatics.tex` | Main manuscript |
| `supplementary_material.tex` | Supplementary PDF (compile separately if needed) |
| `figures/*.png` | All six figures (pre-rendered PNG, 180 dpi) |
| `references.bib` | Bibliography |
| `oup-authoring-template.cls` | OUP document class |
| `oup-abbrvnat.bst` | Citation style |
| `data/` | Small CSV/JSON sources to regenerate figures |

## Figures

All `\includegraphics` paths use explicit `.png` extensions under `figures/`:

- `fig1_error_rate_per_ccd.png` — main Fig. 1a
- `fig2_deposition_year_vs_errors.png` — main Fig. 1b
- `fig6_ramachandran_expanded.png` — main Fig. 2
- `fig3_signed_volume_distribution.png` — supplementary
- `fig4_bland_altman.png` — supplementary
- `fig5_ccd_heatmap.png` — supplementary

Regenerate from bundled data:
```bash
python3 figures/generate_figures.py
```

## Supplementary PDF

In Overleaf, change the main document to `supplementary_material.tex` and recompile, or add it as a second project.

## Page limit

After compiling, check the main PDF page count (Bioinformatics Original Paper limit: **7 pages**).
