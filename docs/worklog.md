# Work log, figures-v2

## 2026-10-03, end of the calibration session
Done: fresh tree; figures/shared/primitives.tex (crescent atoms, rods, solid and ghost molecule modes,
leader labels, zoom bubble, plates, level rules, thick-tube pic); studies/primitives-proof.tex.
Verdict on the proof sheet (build/studies/primitives-proof.png): the molecule primitive passes
(solid methylamine, ghost copy, water, K and Rb discs). The thick-tube pic FAILS and must be rebuilt
before it is used: the top annulus is drawn complete instead of leaving the wedge sector open; the
cut faces are drawn in the wrong order and read as a flat tinted rectangle; the ghost band inside the
wedge is missing; version points fall inside the opening. The working reference for the correct
construction and order is studies/tube-study.tex and studies/harter-idiom.tex (inner wall, outer wall
pieces, ghost band with the tu path, cut faces hatched, top face as annulus minus the sector, then
points and leaders). Rebuild the pic from that order, prove it on the sheet, then proceed to the
prototypes (figures 6, 5, 10).
Nothing in this tree is yet at the standard of section 8. Keep iterating.

## 2026-10-04, continuation after the author's review round
- Rod junctions computed (bonds leave the far sphere on the projected junction curve); digits inside hydrogens.
- Figures 1, 5, 6, 7, 9, 10 reworked after feedback; figures 2, 3, 4, 8, 11, 12 and the guide drawn from
  alternative studies (studies/alternatives, docs/illustration-notes.md); twelve-approach exploration of the
  solid and eight-approach exploration of the molecule shading.
- Reference material kept in the repository (manuscript export, handoff package, plate scans).
- build.py: proof sheet, squint sheet, grayscale proof, deliverable PDFs, collision report; scaffold rebuilt.
- Open for the author: section order of the export against the specification (docs/author-decisions.md).

## 2026-10-04, late: palette anchored, cover accepted, plan v3

Palette pitched on ten cover versions and anchored (see `docs/PLAN-v3.md`); actions moved to ink line styles; the
cover made cohesive (white sheets, lighter hatch, headings only). Figure 7 drafts A, B, C in
`studies/alternatives/round3/`. Study helper used this session (recreate in the scratchpad if missing):

```
#!/bin/bash   # study.sh file.tex [dpi]: compile into build/alt with the shared and data paths, render a png
ROOT=/home/user/Molecular_States; SRC="$1"; DPI="${2:-110}"
export TEXINPUTS="$ROOT/figures/shared//:$ROOT/figures/data//:$ROOT/$(dirname "$SRC")//:"
cd "$ROOT/build/alt" && pdflatex -interaction=nonstopmode -halt-on-error -output-directory "$ROOT/build/alt" "$ROOT/$SRC" > /dev/null \
  && pdftoppm -r "$DPI" -png -singlefile "$ROOT/build/alt/$(basename "$SRC" .tex).pdf" "$ROOT/build/alt/$(basename "$SRC" .tex)"
```
Next: figure 7 (fix the disc size test in `fig07-common.tex`: compare in cm, not pt), then 8, 9.

## 2026-10-04, end of day

Installed from the round-3 studies: figure 7 (sites on the torsion circle, disc-picture grid, the E plane), figure 8
(spin space as hairline bundles between the two parity columns; the dot stacks and the mosaic were rejected as
tables), figure 9 (unrolled packets in three rows, gold only on the overlap and the A1 share), figure 10 (band,
unroll, chart; slice beside q; both normal directions), figure 11 (version pairs side by side, t climbs, b hops, the
two orders end apart), figure 12 (labels only, red marks, gold shear line), figure 4 (the reaction plate with gold and
teal strands), figure 5 (generator line styles over a gold chain halo). `checks/verify_spin.py` added to the build.
Full build, manuscript and zip refreshed at the end (see git log).
