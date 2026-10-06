# Quantum State Spaces of Nonrigid Molecules and Their Complexes

The manuscript `docs/manuscript.tex` and its fourteen plates: the cover guide, twelve numbered plates (one per
section) and the unnumbered plate of the closing example. Each plate is a standalone TikZ document drawn from
computed, independently checked data.

## Build
```
python3 build.py              # data, checks, plates, manuscript (build/manuscript.pdf)
python3 build.py --only 06    # one plate (data and checks still run)
```
A failing check stops the build. Requires TeX Live (pdflatex, bibtex, TikZ, PGFPlots, standalone, titlesec),
Python 3 with SymPy and poppler (pdftoppm, pdftotext). GAP, when present, runs an independent cross-check. The checks
work in exact arithmetic (integers, rationals, the cyclotomic field of the cube roots of unity) wherever the quantity
allows.

## Layout
```
docs/manuscript.tex         the manuscript; docs/refs.bib its bibliography
figures/src/figNN-*.tex     one standalone document per plate (fig00 the guide, fig13 the closing plate)
figures/shared/             figurestyle.sty (palette, type) and primitives.tex (lines, atoms, rods, tube, hatch)
figures/data/               written by compute/make_data.py, never by hand
figures/pdf/, figures/png/  the compiled plates, vector and 600 dpi
compute/                    geometry, groups and the data emitter
checks/                     independent checks of what the plates draw and the text states; the collision report
build.py                    the one entry point
```
Every plate is 15.2 cm wide and is included at `\linewidth`, so all plates print at one scale.
The numbers in the closing example of rigid methane are written by compute/make_data.py and recomputed in
checks/verify_methane.py and checks/verify.g.
