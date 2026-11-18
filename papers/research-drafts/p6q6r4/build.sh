#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
if ! command -v pdflatex >/dev/null 2>&1; then
  echo 'This manuscript requires pdfLaTeX from a TeX distribution with the standard AMS, Latin Modern, microtype, and hyperref packages.' >&2
  exit 1
fi
mkdir -p build_logs
for pass in 1 2 3; do
  echo "Building p6q6r4_manuscript, pass ${pass}"
  pdflatex -interaction=nonstopmode -halt-on-error p6q6r4_manuscript.tex >"build_logs/p6q6r4_manuscript_${pass}.txt"
done
echo 'Built p6q6r4_manuscript.pdf.'
