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

## 2026-10-04, evening

Notation migration (docs/notation-migration.md): the approved symbol table applied role by role across the
manuscript, the plates, the style, the exposed data names and the live docs; checks pass; no undefined macros. Then
the author's plate-by-plate round (docs/illustration-notes.md, round 4): cover with the grey bond interval at the
top left (pitch A of three, tuned), figures 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 reworked to the notes; the parts are
roman-numbered sections, the sections subsections 1-12; a notation table before the references. Full build,
manuscript (20 pages) and zip refreshed at each commit; everything pushed to figures-v2.

## 2026-10-04, night: closing entry

References and bibliography compiled into the manuscript (docs/refs.bib, bibtex step in the README recipe); the
manuscript is self-contained (no legacy motif inputs). The formatting pass (round 5) is applied on all thirteen
plates; rho = 1/2 with Mellor et al. 2019 cited. Full build clean, all checks pass, no text collisions; the proof
sheet, previews, figures/pdf and the zip are current; everything is pushed to figures-v2. Open for the author:
nothing on the list; the next round starts from the author's reading of the draft.

## 2026-10-04, night, round 6 (author's notes on the core and the hatch)

Core sheet half as thick everywhere (tube/bw .04; the zoom strip to match); figure 10's stripped core drawn as a thin
teal tube of that thickness with the seam marks, the small solid without paths, the chart's t arrow bowed clear of the
edge; figure 6's miniature white with red dots on top and pale dots on the hidden edge, the grey images without a
ring; "cladding" is "vibrations" on the plates and in the notes; the cover's slab hatched with continuous lines along
its own slanted edge (a tiled pattern breaks at every tile, so hatch is now drawn by \hatchregion).


## 2026-10-04, night, round 7 (the final automated edit before manual control)

Manuscript content items from the outside review applied (projector, species glossary, compensating rotation,
Mellor attribution, reduced torsional model, agenda dispositions, section 7 order, section 8 bullets and table);
compiles to 20 pages with no errors or unresolved references. Plates: no grey lines; vestigial marks and styles
removed; figure 6's Newman row redrawn as two steps plus a reference; the cover's drop through the centre of mass,
dotted; arrowheads and leaders audited on every plate; two readers' lists applied. Full build clean, all checks
pass (the share data now starts at Delta/d = 0 with the limit shares), no text collisions. Handed to the author
for manual control.

## 2026-10-04, late night, the rigid-formulation box (final pass before manual control)

The closing box of Section 12 inserted from the author's draft after verification against Albert et al.
(arXiv:2403.04572v4); four small technical adjustments recorded in docs/author-decisions.md; the Section 6, 11 and
12 pointers repaired; the bibliography entry names the version consulted. Checks: verify.py, verify_symbolic.py,
the GAP cross-check, verify_spin.py and verify_antiprism.py all pass; the manuscript compiles to 21 pages with no
errors, no over- or underfull boxes and no unresolved references. The prose diff is 86 lines (two files).

## 2026-10-04, late night, the closing page

The rigid-formulation box replaced by a one-page explanation with an unnumbered plate (fig13-rigid-recovery): frozen
shape and one version against the free family and six versions, and the correlation diagram A' -> A1 + E (20 = 12 +
8), A'' -> A2 + E (12 = 4 + 8). The rigid weights are emitted by make_data and checked in verify_spin.py (Frobenius).
Notation and References on their own pages. Minimal captions on all figures; the cover without its caption line.
Build clean, 22 pages.

## 2026-10-04, late night, the closing page on methane

The rigid-recovery page rebuilt on rigid methane for a one-to-one check against Albert et al.: a methane geometry and
glyph in the data layer, the closing-page numbers emitted by make_data, an independent check (verify_methane.py:
24 relabellings as proper rotations, 5 A1 + E + 3 T2, the parity weights, the four isomers of Table II, Schmidt rank
3, monodromy orders 1, 3, 3, 12) registered in build.py, the plate fig13-rigid-recovery redrawn, the page text
rewritten. Build clean, 22 pages. Flagged: the Newton Kabsch routine fails on third-turns.
