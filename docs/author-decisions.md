# Decisions taken, and decisions left to the author

## Taken during the build (author-approved or author-directed)
- Visual identity locked: two ink weights, two-tone gouache shadow on spheres (the black crescent was replaced at the author's request), generatrix hatching on curved walls, hatch
  only for cut material, flat gouache fills, one oblique projection, labels on leaders or digits inside hydrogens.
- Rods start on the junction curve where the bond leaves the sphere (author: bonds must not enter the nuclei).
- Versions and transformed molecules drawn solid with digits inside the hydrogen discs (author chose option C).
- Lattices: the whole in grey, the retained part in black, the chain in blue, the same on both sides (figure 5).
- Figure 6: zoom shows the cut face itself, hidden versions as faded dashed dots, Newman row keyed to the paths.
- Figure 10: the band of figure 6 is shown being cut and unrolled; the slice at q is the cut face of figure 6.
- Figure 4: no envelopes; dashed boundaries, thin grouping lines, the 16-subgroup lattice of S.
- Figures 8, 12, front page: chosen from alternative studies (docs/illustration-notes.md).
- Water atlas: delta in units of l with the -l/4 column at -0.25 (author correction).
- DECISION items of PLAN-v2 resolved as: matrix inset kept (1); table replaced by the correlation diagram (8);
  slice kept in figure 10 but as the figure 6 cut face, reconstruction molecule kept (10); formula panel dropped
  (12); guide as one scene without frames.

## Open, for the author
1. Section order. The export's sections 5 to 9 (Symmetry Groups, Methylamine, Versions and Symmetry Species,
   Localized States, Feasible Configuration Space) differ from the figure flow of the approved specification.
   As agreed, the integration reorganizes the sections to the figure flow in four parts: `docs/manuscript-reorganized.tex`
   is the export with `\part` headings, the twelve sections of the specification, and the plates placed; the author's
   text is moved verbatim (the Methylamine text continues the Symmetry Groups and Versions section with no heading of
   its own; Nuclear Spin Weights is a new section holding figure 8 and a comment where the author's text goes);
   `docs/integration.patch` is the diff from the export. It compiles with figures 1 to 12 in section order. The
   storyline across the parts is in docs/storyline.md.
2. Final captions, from docs/caption-notes.md.
3. Whether to keep the schematic K2Rb2 positions or replace them with a computed complex geometry.
4. The illustrative component functions of figure 2 and rho = 0.3 of figure 12, if real values are preferred.

## 2026-10-04, late: palette anchored on the cover; ink first; line styles as a budget

- The cover fixes the palette and every plate draws from it alone: ink, red (oxygen on the water plates, versions
  on the methylamine plates), blue (nitrogen only), teal (the reference family as a pale tint; rubidium), gold
  (retained regions; potassium), white, grey. Green and violet left the palette.
- Plates read as black-and-white ink drawings; colour only where it speeds comprehension or guides the eye;
  surfaces that carry no meaning stay white; walls of colour are harsh.
- Actions carry no colour: solid t, dashed u or tu, dotted a starred operation. Line styles and hatching are a
  budget per plate like colour; hidden lines are dotted because dashed is spent.
- Leaders start and end on the things they join; labels sit on the object they name.
- The cover's sheets lift the same two loops as the solid; the non-abelian point is left to figure 11.
- A plate's content that is a rule is drawn as the rule acting, not as its outcome (figure 8).

## 2026-10-04, night: symbol clashes resolved document-wide

- The packet width is lambda (was sigma, which is the permutation symbol in 31 of its 33 uses): text lines on the
  Gaussian overlap, figure 9, caption notes.
- The second spin basis vector in the ambient-space example is zeta (was eta, which is the amino coordinate
  everywhere else): one sentence of the text, figure 2's caption note.
- Free letters checked before choosing: lambda, zeta, nu and Delta are unused elsewhere.

## 2026-10-04, night: colour scheme locked

The current hues are final: ink 303438, red D56860, blue 397BA8, teal 4F9C97, gold E0B84F, white, grey as thinned
ink. The grayscale study (L* 21, 49, 57, 60, 77; even and wide respacings; a nitrogen nudge to L* 45) was reviewed
and set aside: hues are used pure, never darkened; fades mark only obscured elements.

## Notation migration (approved instructions, 2026-10-04)

Applied as a bounded update across the manuscript source, the plate sources, the shared style, the exposed data
names and the live documentation. Role by role: configuration matrices X, X0 bold; the signed right-action matrix
P_s, the rotation matrices R(r), R_h, R_z, R_y and the chart matrices A_s, M_s bold; the spatial operators U_s and the
projector P_Gamma sans serif; the chart density j_varphi (was script J); the shape tuple a (was q), the normal
coordinates q (was Q), the normal index alpha with n_q coordinates (was a, d); the normal-frame vectors e_alpha(a)
bold; the umbrella coordinate iota (was eta) and the fixed amino offset b_A (was a); the compensating angle
omega_tu (was alpha); the packet width Delta (was lambda, earlier sigma), same formula c = exp(-d^2/8 Delta^2), no
rescaling; H_amb (was H_ambient); H_+, H_- (were superscripts). The record of what was mapped where is
docs/notation-migration.md.

## Round of plate notes (author, 2026-10-04, afternoon)

- Cover: the bond interval [H, B] joins the cover in grey at the top left (pitch A of three, with the lattice smaller
  and a step down and right, a vertex named G, the tube, sheets and their headings slid right by .45 for balance);
  the free-domain label V is dropped everywhere (the quotient is named on the projection; figure 11's base is F/G).
- Structure: the four parts become roman-numbered sections, the twelve sections subsections numbered 1-12 straight
  through; the cover maps the four sections. A notation table stands before the references.
- Actions and spins in ink: figure 3's spin arrows are ink (the coloured ones were hard to see).
- Figure 6 reads as three levels of zoom (C with twenty images, one an exact miniature; F; the cut face), figure 10
  as a sequence (solid, core sheet, chart) with the two fibres over a point drawn as the reconstruction.
- Figures 7, 8 and 9 take less room; figure 8's bundles have no end strokes; figure 9's share chart spans the rows
  with the three ratios written at their dots.

