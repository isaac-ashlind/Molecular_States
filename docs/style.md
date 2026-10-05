# Style of the suite

Reference: the ink plates of Harter (1993); for line and flat colour also Winsor McCay (clean contour, restrained
flat fills, air around every figure). The drawing shows the structure itself; text is bookkeeping. The plates read
as black-and-white ink drawings; colour is spent only where it speeds comprehension or guides the eye. Surfaces that
carry no meaning stay white. Everything below lives in `figures/shared/figurestyle.sty` (palette, `plate`) and
`figures/shared/primitives.tex` (line styles, atoms, rods, tube, band, hatch, zoom bubble).

## Line
- Ink `#303438` for every line. No grey lines, except for greyed-out objects (the cover's lattice, the inactive parts
  of the lattices in figures 4 and 5, faded shapes).
- Weights: `heavy` 0.9 pt (silhouettes, retained structure, action paths); `fine` 0.4 pt (construction, detail);
  0.3 pt (leaders, axes behind a molecule, the cube on the closing plate); hatch 0.22 pt.
- Actions carry no colour; the line style says which, at action weight: solid = t, dashed = u or tu, beaded dots =
  a starred operation (b, E*). Stealth heads only where a path has a direction.
- A hidden part of a contour is the same stroke dashed (`back`: on 1.5 pt, off 1.1 pt), on every plate.
- A hidden or faded object's outline is dashed (`ghost`); a hidden point (every edge at it hidden) is pale in a ghost
  outline, drawn under the front lines. One exception: the undisplaced ghost under the displaced molecule of figure 10
  is a simple faded outline, solid, with the suite's rods and white atoms (author).
- A projection or lift between levels is thin and dotted (`projection`), as on the cover and in figures 11 and 13.
- Leaders (`leader`, 0.3 pt) start on the object and stop short of the label; they cross no other element; where a
  short leader would cross clutter it is made longer instead. No white halos.

## Colour
Gouache, anchored on the cover plate and used exclusively: red `#D56860` (oxygen on the water plates, version points
on the methylamine plates and the rotations of T on the closing plate), blue `#397BA8` (nitrogen only), teal
`#4F9C97` (the reference family as a pale surface tint, the orientation sphere; rubidium at full strength), gold
`#E0B84F` (retained regions: the chart patch, the cell, the base of a covering; potassium), white (cladding,
hydrogen). No other hue appears anywhere in the suite.
- Spheres: two tones of one gouache, the darker as the shadow lower right (light from the upper left), ink outline.
  Version points are such spheres at radius 0.08 cm.
- Rods: white band with ink edges, the shadow edge heavier, starting on the junction curve of the far sphere.
- Curved walls: generatrix lines crowding toward the silhouettes.
- Cut material is hatched in ink at 0.22 pt with continuous lines (45 degrees; the cover's slab follows its receding
  edges). Spacing and phase are set so that an edge parallel to the hatch falls midway between two lines, never on
  one (`\hatchfit`).

## Projection
- Methylamine: one oblique camera (azimuth 50, elevation 25). Solids: ellipse ratio 0.42. The K sphere of figure 12
  at the higher elevation 0.6, so the ring planes open.
- The closing plate shares one right-handed camera between methane and the orientation ball, so bond 1 points along
  the arrow. (The suite's `G.camera` builds a left-handed frame; see `HANDOFF.md` on why the methylamine plates are
  unaffected in what they claim.)

## Type
- Computer Modern. Labels 9 pt (`lab`), notes 8 pt (`sub`); the guide's part headings 9.5 pt bold. Smaller sizes only
  inside dense marks: digits in hydrogen discs 6.5 pt; figure 5's cover letters 6.5 pt and lattice labels 7.5 pt;
  figure 11's graph nodes and figure 1's cells 7 pt.
- Every plate is 15.2 cm wide and is included at `\linewidth`, so all plates print at one scale: labels at 9.8 pt,
  notes at 8.7 pt.
- A label sits on what it names, with air, never on a line, an edge, a hatch or another label.
- Text under an object is centred on it; a row or corner title is left-aligned at the margin.
- In-plate text is lowercase phrases without full stops. Captions in the manuscript are minimal.
- Headings in the manuscript: a period after the number, none after the title.

## Depth
2D on purpose where the idea is 2D (chart, lattice, correlation diagram, radial plots); depth only where it carries
intuition (molecules, the solid, the sphere ladder, the orientation ball, the sheets).
