#!/usr/bin/env bash
# Build Bioinformatics submission PDFs and upload zip.
set -euo pipefail
cd "$(dirname "$0")"

MAIN=chiralfold_bioinformatics
SUPP=supplementary_material
ZIP=ChiralFold_Bioinformatics_Submission.zip

echo "==> Building main manuscript..."
pdflatex -interaction=nonstopmode "$MAIN.tex" >/dev/null
bibtex "$MAIN" >/dev/null || true
pdflatex -interaction=nonstopmode "$MAIN.tex" >/dev/null
pdflatex -interaction=nonstopmode "$MAIN.tex" >/dev/null

echo "==> Building supplementary material..."
pdflatex -interaction=nonstopmode "$SUPP.tex" >/dev/null
pdflatex -interaction=nonstopmode "$SUPP.tex" >/dev/null

echo "==> Creating submission zip..."
rm -f "$ZIP"
zip -r "$ZIP" \
  "$MAIN.tex" "$MAIN.pdf" \
  "$SUPP.tex" "$SUPP.pdf" \
  references.bib \
  oup-authoring-template.cls oup-abbrvnat.bst \
  COVER_LETTER.txt HIGHLIGHTS.md README_SUBMISSION.md \
  figures/

echo "==> Done."
echo "Main PDF: $(pwd)/$MAIN.pdf ($(pdfinfo "$MAIN.pdf" 2>/dev/null | awk '/Pages/ {print $2}') pages)"
echo "Supplementary: $(pwd)/$SUPP.pdf"
echo "Upload zip: $(pwd)/$ZIP"
