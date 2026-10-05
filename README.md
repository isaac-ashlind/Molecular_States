# Quantum State Spaces of Nonrigid Molecules and their Complexes

The manuscript `docs/manuscript.tex` and its fourteen plates: the cover guide, twelve numbered plates (one per
section) and the unnumbered plate of the closing example. Each plate is a standalone TikZ document drawn from
computed, independently checked data.

## Build
```
python3 build.py              # data, checks, plates, manuscript (build/manuscript.pdf)
python3 build.py --only 06    # one plate (data and checks still run)
```
A failing check stops the build. Requires TeX Live (pdflatex, bibtex, TikZ, PGFPlots, standalone, titlesec),
Python 3 and poppler (pdftoppm, pdftotext). numpy, SymPy and GAP run further checks; the build says when it skips one.

## Layout
```
docs/manuscript.tex         the manuscript; docs/refs.bib its bibliography
figures/src/figNN-*.tex     one standalone document per plate (fig00 the guide, fig13 the closing plate)
figures/shared/             figurestyle.sty (palette, type) and primitives.tex (lines, atoms, rods, tube, hatch)
figures/data/               written by compute/make_data.py, never by hand
figures/pdf/, figures/png/  the compiled plates, vector and 600 dpi
compute/                    geometry, groups and the data emitter
checks/                     independent checks of every number the plates draw; the text-collision report
build.py                    the one entry point
```
Every plate is 15.2 cm wide and is included at `\linewidth`, so all plates print at one scale.
