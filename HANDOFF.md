# Hand-off

The manuscript and its plates are complete, built from one command, and checked. This note says what is where, what
is verified, how to change things, and what is left to the author.

## State
- `docs/manuscript.tex` builds to 24 pages (`build/manuscript.pdf`): the cover guide, the twelve sections in four
  parts with one plate each, the three closing pages on rigid methane, the notation and the references. No undefined
  references, no over- or underfull boxes.
- Fourteen plates in `figures/src/`, compiled to `figures/pdf/` (vector) and `figures/png/` (600 dpi). Every plate is
  15.2 cm wide and is included at `\linewidth`, so all print at one scale (labels 9.8 pt, notes 8.7 pt).
- Branch `figures-v2`. The backup taken before the final sweeps is commit `7765ac7` (and `ba92480` after the heading
  change); the source zip and PDF of both were delivered. The tag `v2-backup-before-sweeps` exists locally only: the
  remote refuses tag pushes from this environment.

## Build
```
python3 build.py              # data -> checks -> plates -> previews -> numbering check -> proofs -> manuscript
python3 build.py --only 06    # one plate
```
Outputs in `build/`: `manuscript.pdf`, `proofsheet.pdf` (every plate at print width, then the squint sheet),
`proofsheet-gray.pdf`, `previews/` (150 dpi colour and grayscale), `collisions.txt`; the plates themselves go to
`figures/pdf/` and, at 600 dpi, `figures/png/`. A failing check stops the build.

## Where to change what
- **Text**: `docs/manuscript.tex`. Banner comments mark the four parts, the twelve sections and the closing pages.
  The closing pages have explicit page breaks (before item 5 and before Table A); if their text grows, move those.
- **A plate**: `figures/src/figNN-*.tex`, then `python3 build.py --only NN`. The visual vocabulary (line styles, atoms,
  rods, tube, hatch, zoom bubble) is in `figures/shared/primitives.tex`, the palette in `figurestyle.sty`; the rules
  are written out in `docs/style.md`. What each plate shows and what it declares schematic is in
  `docs/caption-notes.md`.
- **A number or a coordinate**: `compute/make_data.py` emits everything the plates draw into `figures/data/`; the
  checks in `checks/` recompute it independently and compare. Never edit `figures/data/` by hand.

## What is verified
`docs/verification.md` has the full list. In short: group orders, cosets, the subgroup intervals and covers, spin
weights and species (stdlib, SymPy and GAP independently); the band's symmetry; every number and construction on the
closing pages, from explicit matrices (`checks/verify_methane.py`); text overlaps on every plate. All pass.

## Rules the author set (the short list; `docs/style.md` and `docs/author-decisions.md` have the rest)
- Ink lines only, no grey lines except greyed-out objects; colour only where it carries meaning.
- Actions by line style at action weight: solid t, dashed u or tu, beaded starred. Hidden contours dashed in their own
  stroke; projections dotted; leaders thin, starting on the object, stopping short of the label, crossing nothing.
- Labels on what they name with air; text under an object centred; in-plate text lowercase without full stops.
- Minimal captions; headings with a period after the number only.

## Open items for the author
1. **Camera handedness of the methylamine glyphs.** `compute/geometry.py`'s `camera()` builds a left-handed frame, so
   every glyph drawn through it is a mirror image: all 3D methylamine glyphs (cover, figures 5, 10, 11). Nothing the
   reader is told is affected: the mirror is uniform, the inversion E* commutes with every permutation and rotation,
   so every stated relation (bX0 = R_b X0 with R_b = R_y(pi), the lifts to tX and bX, the parity of the normal
   displacements under tu) holds the same way in the mirror copy; no mirrored glyph sits beside a true view of the
   same configuration (figure 6's Newman projections are true views and figure 6 has no 3D glyph); and no turning
   sense is stated for a mirrored glyph. The closing plate uses its own right-handed camera. To make every glyph a
   true view, change `out = cross(up, right)` to `cross(right, up)` in `camera()` and re-check each plate's look
   (depth order flips; layouts were tuned by eye).
2. **`geometry.kabsch_rotation`** (a Newton polar iteration) does not converge for third-turns; it is exact for the
   small turns and half-turns the methylamine plates use. The methane check fits by SVD instead.
3. **Captions** are minimal by decision; the disclosures a careful caption may want (schematic elements, declared
   approximations) are listed per plate in `docs/caption-notes.md`.
4. **`reference/`** holds the author's materials and was not touched.
5. **Two tight spots left by the last label pass** (every leader and label was measured; these need geometry, not a
   nudge): in figure 1 the 7 pt entries of the column strips sit 1.3 to 2 pt from the cell rules (wider cells or
   smaller text would fix it); in figure 11 the labels of the bookkeeping graph nearly touch their circles (bigger
   nodes would fix it). Two leaders cross more than one line because their dot sits inside shading or hatch:
   figure 6 t²uH (five thin wall lines, the fewest possible) and sF (the hatched slab).

## Last changes (2026-10-05)
Each is logged in `docs/author-decisions.md`; every change is its own commit on `figures-v2`, so each diff can be read
alone.
- Before the deep review: figure 4's strands one weight and its arrows centred; figure 5's fragments turn with their
  clips; figure 10's ghost a simple outline; the closing plate's hidden points pale; a label pass on eight plates.
- Notation (1a65f72, 5d0fb69, 9f64980): the closing pages write the vibrations as q in R^9, not the disc D^9 (which
  clashed with the Wigner D^J); a full audit then gave every symbol one meaning (omega only an angle, packets eta,
  loops lambda, s' in the action rule, the proton spin space (C^2)^{x5}, Section 2's components f_xi and f_zeta, Table
  A's classes by relabellings), defined the rest at first use, and completed the notation page.
- Plates (244892e to 24fd9a7): a seven-inspector, seven-skeptic audit; 46 confirmed findings applied plate by plate.
