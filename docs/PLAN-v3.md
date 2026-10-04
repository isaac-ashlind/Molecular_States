# Plan v3: the road ahead (written 2026-10-04, robust to interruption)

Branch `figures-v2`. Everything below is committed; resume from this file. The author is the physicist, the
assistant is the scientific illustrator. No swarms: one coordinator, at most two bounded subagents (reviewers).

## Anchored decisions (do not reopen)

- **Palette, exclusive and locked** (`figures/shared/figurestyle.sty`): ink `#303438`; red `#D56860` (oxygen on the water
  plates, version points on the methylamine plates); blue `#397BA8` (nitrogen only); teal `#4F9C97` (the reference
  family as a pale surface tint; rubidium at full strength); gold `#E0B84F` (retained regions; potassium); white
  (cladding, sheets, hydrogen); grey as thinned ink. Any palette colour not already spent on a plate is available
  there, if it does not conflict with the storyline.
- **Ink first.** Plates read as black-and-white ink drawings; colour only where it speeds comprehension or guides
  the eye. Surfaces that carry no meaning stay white (the sheets over V are white).
- **Actions carry no colour.** Solid = t, dashed = u or tu, dotted (beaded) = a starred operation (b, E*).
  Styles `tpath`, `upath`, `bpath`, `tlift`, `ulift`, `blift` in `figures/shared/primitives.tex`.
- **Cover** (`figures/src/fig00-guide.tex`): accepted composition. Slab C with lighter hatch, molecule as a point,
  solid F lifted out with the teal core sheet, red versions, gold packet, t solid and tu dashed; V with the two
  loops; four white sheets (H, uH, tH, tuH) with the same two loops lifted, ending apart; four headings with
  section ranges only. Monodromy is shown by lifts that do not close; the non-abelian point is left to figure 11.
- Versions drawn solid with digits inside hydrogen discs; two-tone gouache sphere shading; bonds start on sphere
  surfaces; no labels on hatching; a zoom has the same shape as what it points at; hatch marks cut material only.
- Sections in four parts in figure order, no "Methylamine" heading (`docs/manuscript-reorganized.tex`).

## The loop on every plate

1. One sentence: what the reader must see.
2. Three distinct approaches as thumbnails, side by side, in the palette (`studies/alternatives/round3/`).
3. Pick or merge; refine at plate size; crops at 220 dpi; `checks/collisions.py`; grey proof.
4. Word budget: labels plus at most one short line; explanation goes to `docs/caption-notes.md`.
5. Adopted elements propagate (version dots, torsion circle, band and chart, disc pictures, line grammar).
6. Record the critique in `docs/illustration-notes.md`; commit after every plate.

## Status (end of 2026-10-04)

Done this round: cover accepted; figures 7, 8, 9, 10, 11, 12, 4 and 5 redrawn and installed (studies and critiques in
`studies/alternatives/round3/` and `docs/illustration-notes.md`); palette, ink-first rule and line-style budget
written into the style files. Open: the author's review of the installed plates; touch-ups to figures 1 to 3;
the reviewer-subagent close-up pass; the final QA list below.

## Per plate (order: 7, 8, 9; then 10, 11, 12; then 4, 5; then 1 to 3 touch-ups)

- **7 species.** Sites on the torsion circle as in 9 and 12; t the arc, b the mirror line through H; modes as disc
  pictures on the sites; grid of each mode with its t image and b image (A1 unchanged, E mixes). Alternative: the E
  plane with the three images of v2 at a third turn and the mirror axis. No matrices. Drafts A, B, C exist.
- **8 spin weights.** The 32 spin states as one strip partitioned 12, 4, 16, joined to the three spatial level
  rules by solid (even) and dashed (odd) lines; weights read as lengths. Alternative: dot columns above levels.
- **9 localized states.** Unrolled tau line, three packets with the overlap shaded, rows for narrow, drawn, broad,
  each with its A1 share as a bar; chart small with the rows marked. Alternative: radial profiles as in 12.
- **10 position representation.** Unrolling as a sequence (band, cut, flat), chart with fine cells, reference cell
  gold, six red points, two ink paths, formulas off the plate, slice zoom and one molecule with room.
- **11 monodromy.** Twelve white sheets in six pairs; t then b and b then t ending apart; endpoints as gold rings,
  never a dot over text; the six-sheet graph with end nodes gold; X and bX magnifier.
- **12 momentum.** Sphere ladder, two lattices, radial modes; the torsion circle identical to 7 and 9 with red
  version marks; shear line gold; paragraphs cut to labels and level rules.
- **4 reaction.** G_in left; gold and teal strands continuous through S; reaction-path panel at right in the same
  colours; text block removed.
- **5 groups and versions.** Gold chain; covers by generator line style; mirror glyph dropped; versions row given room.

## Acceptance before "done"

`python3 build.py` clean; 12 figures plus the guide; `build/collisions.txt` empty of overlaps; grey proof legible;
`docs/manuscript-reorganized.tex` compiles with the guide on page one; `dist/molecular-figures-v2.zip` refreshed;
docs updated (style, caption notes, illustration notes, verification, author decisions, worklog); pushed.

## How to resume

```
python3 build.py                                  # data, checks, figures, previews, proof sheets
scratch/study.sh <file.tex> [dpi]                 # see docs/worklog.md for the one-line helper
python3 checks/collisions.py build/figures/fig07-symmetry-species.pdf
```
Studies and their critiques: `studies/alternatives/round3/`, `docs/illustration-notes.md`.
