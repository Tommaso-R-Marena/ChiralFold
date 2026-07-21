# ChiralFold Bioinformatics Overleaf package

Use `ChiralFold_Bioinformatics_Overleaf.zip` for one-click Overleaf upload.

1. Download the `Overleaf/` folder or the zip file from this repository.
2. In Overleaf, choose **New Project** → **Upload Project**.
3. Upload `ChiralFold_Bioinformatics_Overleaf.zip`.
4. Set the main document to `chiralfold_bioinformatics.tex` and recompile with pdfLaTeX + BibTeX.

Optional local previews (already compiled under TeX Live) are in
`compiled_preview/`:

- `chiralfold_bioinformatics.pdf` — main manuscript (~6 pages; Bioinformatics Original Paper limit is 7)
- `supplementary_material.pdf` — supplementary PDF

Also in this folder: `COVER_LETTER.txt`, `HIGHLIGHTS.md` (not required for Overleaf compile).

To refresh the zip from source, run `scripts/build_overleaf_package.sh` from the repository root.
