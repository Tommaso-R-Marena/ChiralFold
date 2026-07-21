#!/usr/bin/env bash
# Fail if LaTeX body paragraphs spill badly off the page (>5pt overfull).
# Decorative OUP header rules (~262pt) are ignored.
set -euo pipefail
cd "$(dirname "$0")"
THRESHOLD="${1:-5}"

check_log() {
  local log="$1" name="$2"
  local bad=""
  while IFS= read -r line; do
    [[ "$line" == *"261.76535"* ]] && continue
    if [[ "$line" =~ \(([0-9.]+)pt\ too\ wide\) ]]; then
      local w="${BASH_REMATCH[1]}"
      if awk -v w="$w" -v t="$THRESHOLD" 'BEGIN { exit (w+0 > t+0) ? 0 : 1 }'; then
        bad+="${line}"$'\n'
      fi
    fi
  done < <(grep 'Overfull \\hbox' "$log" || true)
  if [[ -n "$bad" ]]; then
    echo "==> $name: overfull paragraphs exceeding ${THRESHOLD}pt:"
    printf '%s' "$bad"
    return 1
  fi
  echo "==> $name: OK (no body overfulls > ${THRESHOLD}pt)"
}

pdflatex -interaction=nonstopmode chiralfold_bioinformatics.tex >/dev/null || true
bibtex chiralfold_bioinformatics >/dev/null 2>&1 || true
pdflatex -interaction=nonstopmode chiralfold_bioinformatics.tex >/dev/null || true
pdflatex -interaction=nonstopmode chiralfold_bioinformatics.tex >/dev/null || true
pdflatex -interaction=nonstopmode supplementary_material.tex >/dev/null || true
pdflatex -interaction=nonstopmode supplementary_material.tex >/dev/null || true

check_log chiralfold_bioinformatics.log "main"
check_log supplementary_material.log "supplementary"
