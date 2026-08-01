#!/usr/bin/env bash
# Build the self-contained Bioinformatics Overleaf upload package at repo root.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SUBMISSION_DIR="$ROOT/paper/submission/bioinformatics"
OUT_DIR="$ROOT/Overleaf"
ZIP="ChiralFold_Bioinformatics_Overleaf.zip"

"$SUBMISSION_DIR/package_overleaf.sh"
mkdir -p "$OUT_DIR"
cp "$SUBMISSION_DIR/$ZIP" "$OUT_DIR/$ZIP"

# Refresh the committed local previews when a TeX toolchain is available, so
# Overleaf/compiled_preview/ never drifts from the zip it is meant to preview.
if command -v pdflatex >/dev/null 2>&1; then
  echo "==> Compiling local previews..."
  (
    cd "$SUBMISSION_DIR"
    pdflatex -interaction=nonstopmode chiralfold_bioinformatics.tex >/dev/null 2>&1 || true
    bibtex chiralfold_bioinformatics >/dev/null 2>&1 || true
    pdflatex -interaction=nonstopmode chiralfold_bioinformatics.tex >/dev/null 2>&1 || true
    pdflatex -interaction=nonstopmode chiralfold_bioinformatics.tex >/dev/null 2>&1 || true
    pdflatex -interaction=nonstopmode supplementary_material.tex >/dev/null 2>&1 || true
    pdflatex -interaction=nonstopmode supplementary_material.tex >/dev/null 2>&1 || true
  )
  mkdir -p "$OUT_DIR/compiled_preview"
  for pdf in chiralfold_bioinformatics supplementary_material; do
    if [ -f "$SUBMISSION_DIR/$pdf.pdf" ]; then
      cp "$SUBMISSION_DIR/$pdf.pdf" "$OUT_DIR/compiled_preview/$pdf.pdf"
    fi
  done
else
  echo "==> pdflatex not found; leaving Overleaf/compiled_preview/ unchanged"
fi

# Keep the standalone cover letter / highlights copies in step with the source.
for f in COVER_LETTER.txt HIGHLIGHTS.md; do
  [ -f "$SUBMISSION_DIR/$f" ] && cp "$SUBMISSION_DIR/$f" "$OUT_DIR/$f"
done

echo "Wrote $OUT_DIR/$ZIP"
