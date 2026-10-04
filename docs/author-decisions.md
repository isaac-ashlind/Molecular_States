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

## Second review round and the frame parameter (2026-10-04, evening)

- Two bounded reviewers (one per half of the suite) listed small formatting defects; all concrete items were applied
  (docs/illustration-notes.md, round 5). Declined: recolouring figure 2's components (teal and gold are the approved
  dual pair there) and the edge-on cut face of the small solid in figure 10 (geometry, not a defect).
- Conventions fixed by that round: version points radius .08 on every plate; chart references (loci, zero lines,
  norm curves) fine solid lines, never a dashed or dotted style (grey then, ink since the round below); a leader
  starts on the edge of the point it joins;
  an arrow shaft starts on an atom's silhouette; the same molecule at one scale within a plate.
- The frame-twist parameter is rho = 1/2, the symmetric (internal-axis) frame of Mellor, Yurchenko, Mant and
  Jensen, Symmetry 11, 862 (2019), cited in the text; figure 12's twisted labels are the half-integers m + K/2.
  Numbers on the plates and in the text are written as exact fractions and roots, never as decimals.
- Plate sizes: a plate shares its page with the prose (at most about two thirds of the text height) or takes the
  page, which only figure 5 does; floats are no longer barred at subsections, so the prose flows under them.


## The final automated edit (2026-10-04, night; the author's go, then manual control)

- Manuscript content from the outside review, applied as drafted: the projector written with chi_Gamma(g^-1); the
  notation table restricts the species to G6; the compensating-rotation convention stated once (r A_s(a) is the
  orientation whose matrix is R(r) A_s(a)); rho = 1/2 attributed as the analogue for one methyl top of Mellor,
  Yurchenko, Mant and Jensen's symmetric frame; section 12's product model stated as a reduced torsional model with
  iota held at its reference value; the agenda bullets of sections 5, 6 and 12 say where their questions are settled;
  section 7 reordered (species, projector, induction, C[G6/H], inclusions and double cosets); section 8 written out
  with the decomposition 12 A1 + 4 A2 + 8 E, the pairing rule and the weights table.
- No grey lines (author): every line is ink, varied by weight and style; grey stays only on greyed-out objects (the
  cover's lattice, the inactive parts of the lattices in figures 4 and 5, the undisplaced shape in figure 10, hidden
  version points). Converted: the loci of figure 1, the zero line and density curve of figure 2, the E share of
  figure 9, the projections of figure 11 (now dotted like the cover's), the zero line and mode circles of figure 12,
  the inner-wall lines of every tube. Kept grey, as a judgement to be confirmed: the hatch of cut material (a
  texture, lightened in the cover round the author saw).
- Vestigial things removed: figure 2's three unlabelled dots on the density curve (they marked the samples' values,
  which the shared theta scale already gives); the version dots on figure 6's miniature (author); the unused v1
  styles in figurestyle.sty and primitives.tex (contour, guide, operation, inclusion, panel, note, frame, hatch,
  tag, region), an unused cell style in figure 1, an unused dashed style in figure 8, a shadowed definition in
  figure 6; stale header comments (lambda, sigma for the packet width).
- Figure 6's Newman row (author): t is one relabelling, H straight to tH, so the row shows the two steps as one
  arrow each from H (tu to tuH, t to tH) and the eclipsed shape apart at the right as a reference, the midpoint in
  tau of the path of t, not a version.
- The cover's drop from the molecule to its point of C runs through the centre of mass (the glyph is mass-centred,
  so the vertical through the pic origin) and is dotted like the lift (author).
- Arrowheads: a path's head stops about 3pt short of the dot it points at (cover, figures 6 and 10, the chart paths),
  never under it; a head meets the silhouette it points at (figure 1's right P arrow meets the top hydrogen of
  RXP, figure 4's in-arrows end on the complex's edge like the out-arrows); labels sit clear of the barbs.
- Two close readers' lists (about ninety items, scratchpad record readers-round3.md) applied with judgement. Declined:
  one scale for the water glyphs of figure 1 (the atlas must be small, the panel large); raising the cover's lattice
  and heading III by 8 mm (the author set that balance); figure 5's version row at a smaller scale (one scale per
  plate); labelling figure 9's chart dots instead of the axis (the ratios are now ticks).

## The rigid-formulation box (2026-10-04, night; the author's final prompt before manual control)

- One box closes Section 12, before the notation table: the author's approved draft, verified against Albert,
  Kubischta, Lemeshko and Liu, arXiv:2403.04572v4 (their theorem: the rotational state space of a G-symmetric
  isomer is the induced representation Gamma_rot up SO(3); their condition (1): Gamma_rot x Gamma_nuc contains the
  spin-statistics irrep; their assumption of decoupled, trivially transforming electronic and vibrational factors;
  position states over SO(3)/G with the fiber transforming by Gamma_rot(g^-1)).
- Technical adjustments to the draft, each the smallest that correctness needed: (i) H_rot is named as the group
  of the rotations R_h of Section 5 for unstarred h (Albert's proper rotational symmetry group); (ii) the transport
  of the spin species and chi_stat to H_rot is through h -> R_h^-1, because with the right action s.X = X P_s the
  map h -> R_h reverses products (R_{h1 h2} = R_{h2} R_{h1}); the other reading conjugates the characters, so the
  two agree for real characters, which is every character in this paper; (iii) "the induced rotational measure" is
  glossed as the Haar measure dr of Section 10 on the orbit; (iv) one sentence restricts the description to a
  nonlinear reference and points linear configurations to the axial-redundancy remark of Section 10.
- The bibliography entry for Albert et al. now names the version consulted (arXiv:2403.04572v4, 13 January 2026),
  read from the document's own stamp; nothing else was inferred.
- Section 6 no longer claims that shrinking F to a neighbourhood of SO(3)X0 recovers the rigid limit; it points
  to the box. Section 11's unfinished instruction is a pointer to the box. Section 12's opening names its example
  as a reduced torsional model, not a transform of the full torsion-umbrella model.
- Placement: the box cannot share page 18 with figure 12, so it opens page 19 (a \clearpage before it, with a
  measured \Needspace guard that keeps it whole wherever it lands later); page 18 holds figure 12 with blank space
  beneath. Left for manual work: that blank, and any other float interruption; no layout pass was made.
- Not changed, for the author: Section 12's reduced torsional model stays; Sections 7 and 8 untouched; figures
  untouched. Observation only: for methylamine H = {E, b} with b starred, so H_rot is trivial (an asymmetric top in
  Albert's classification) and the whole content of H in the rigid limit is the parity condition carried by b.

## Last plate notes before manual control (2026-10-04, late night)

- Figure 12's sphere: the vector J starts at a marked centre, ends on the front of the K = 1 ring, and the ring's
  radius from the axis to the tip is dotted so the z component reads as K = 1; the sphere is seen from a higher
  elevation than the tubes (E = .6) so the ring planes open. A tilt of the body axis was tried and rejected by the
  author as needless complexity; the axis stays upright.
- Figure 10's seam is a slit (two cut edges, white between) through the back wall with cut marks; H is named outside
  the sheet on a leader. Leaders carry no white halo (author); where a label would sit in clutter the leader is made
  longer instead.
- The hatch of cut material is ink, thinner (.22pt), on the cover and on every plate that uses it.
- The appendix question (a rigid-recovery appendix with a figure) is answered in the hand-off note and not started.
