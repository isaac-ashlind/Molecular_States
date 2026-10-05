#!/bin/bash
# compile docs/manuscript.tex against figures/pdf (pdflatex, bibtex, pdflatex x2) into build/ms; copy to build/manuscript.pdf
set -e
ROOT=/home/user/Molecular_States; cd "$ROOT"; mkdir -p build/ms
for pass in 1 b 2 3; do
  if [ "$pass" = b ]; then (cd build/ms && BIBINPUTS="$ROOT/docs//:" bibtex manuscript > bibtex.stdout 2>&1) || { tail -5 build/ms/bibtex.stdout; exit 1; }
  else pdflatex -interaction=nonstopmode -halt-on-error -output-directory build/ms docs/manuscript.tex > build/ms/ms.stdout 2>&1 || { grep -n "^!" -A3 build/ms/manuscript.log | head; exit 1; }; fi
done
cp build/ms/manuscript.pdf build/manuscript.pdf
pdfinfo build/manuscript.pdf | grep Pages
echo "bad boxes: $(grep -c 'Overfull\|Underfull' build/ms/manuscript.log)  undefined refs: $(grep -c 'undefined' build/ms/manuscript.log)"
