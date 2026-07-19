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

echo "Wrote $OUT_DIR/$ZIP"
