#!/usr/bin/env bash
# Create a self-contained zip for Overleaf upload (compile without extra setup).
set -euo pipefail
cd "$(dirname "$0")"

ZIP=ChiralFold_Bioinformatics_Overleaf.zip

echo "==> Regenerating figures (optional; uses bundled data/)..."
if python3 figures/generate_figures.py 2>/dev/null; then
  echo "    Figures OK"
else
  echo "    Skipped figure regen (using committed PNGs)"
fi

echo "==> Creating Overleaf zip..."
rm -f "$ZIP"
zip -r "$ZIP" \
  chiralfold_bioinformatics.tex \
  supplementary_material.tex \
  references.bib \
  oup-authoring-template.cls \
  oup-abbrvnat.bst \
  figures/*.png \
  figures/generate_figures.py \
  data/ccd_code_coverage_summary.csv \
  data/error_table_verified.csv \
  data/molprobity_comparison.json \
  data/ramachandran_279struct_chainfix_comparison.csv \
  README_SUBMISSION.md \
  OVERLEAF.md

echo "==> Done: $(pwd)/$ZIP"
echo "Upload this zip to Overleaf (New Project -> Upload Project)."
echo "Set main document to: chiralfold_bioinformatics.tex"
echo "Compile with pdfLaTeX + BibTeX (Overleaf default)."
