# Plan v2: a fresh figure repository for the manuscript, built for Overleaf

Status: proposal for the author (2026-10-03). Everything before this document is calibration.
Carried forward unchanged: `compute/` (geometry, groups, data emission), `checks/` (stdlib, sympy,
GAP, CVD), the generated data, and the identity fixed from the author's Harter plates. Every figure
source is rewritten from a blank page.

## 1. Governing idea: carving from the top down
Harter builds exquisite models from the ground up. This manuscript carves the same shapes from the
top down: the full configuration space, the ambient Hilbert space and the complete group come first;
constraints carve the physical space, the feasible region, the species and the localized states out
of them; the last three figures are views of the carved object. The suite reads in that order and
every plate shows the whole before the part.

Visual grammar of carving, one meaning per device, suite-wide:
- the whole, drawn first and lightly, as an outer contour (C, H_ambient, S);
- the retained part in flat gouache with a heavy contour (H in H_ambient, F in C, G in S, G/H in S/H);
- what a constraint removes drawn as cut material, regular 45-degree hatch (hatch never means shading);
- nested outlines read from the outside in (Harter fig. 4.1.5, inverted in reading direction);
- ghosts (dashed, uncoloured) for rotated, relabelled or removed copies.

## 2. Identity (locked)
Ink: two weights only, 0.9 pt silhouettes, principal arrows and section edges; 0.4 pt construction,
leaders, hidden lines (dashed) and hatch (0.3 pt). No grey linework. Shading: solid black crescent on
spheres, solid black band bent around tubes. Gouache: opaque flat fills, ink #303438 (also carbon),
H white, O #D56860, N #397BA8, K #91A653, Rb #AD8BC0, gold #E0B84F for marked points and the
distinguished region; coloured lines only for the t path (blue, solid) and the tu or b path (red,
dashed). 3D: one oblique projection for the whole suite (eye elevation 22 degrees, ellipse ratio
0.37), light from the upper left, hidden lines dashed, rods as double lines occluded by spheres,
molecules at 0.95 cm per angstrom with ink labels on straight leaders. 2D wherever geometry is not the
point. Type: 9 pt serif math, horizontal, 1 mm clear, never over a line; panel letters plus at most
three words; one sub-line of at most twelve words; all other words are the author's captions.
Composition: one idea per panel, at most three panels, object above and bookkeeping beneath joined
by fine leaders, magnifying circles for local enlargements, dotted named loci for special sets.
Citation: a composition adapted from a Harter plate is disclosed in the caption notes as
"after Harter (1993), Fig. x.y.z"; the BibTeX entry is supplied for the author.

## 3. The twelve plates and the guide
Each entry: the one idea; the drawing; what is computed; what is disclosed. Panel contracts are the
approved ones; items marked DECISION would change a contract and wait for the author.

1. Configuration space. Idea: a configuration is a centred matrix; rotations act on the left,
   relabellings on the right, and they commute. Drawing: the water atlas as the carved family
   (25 glyphs from one parameter table, mass-centre cross, dotted equal-bond locus and linear boundary,
   their intersection circled), with the commuting square beside it and the matrix as a small inset.
   Computed: all glyphs, the matrix, X m = 0, P_st = P_t P_s. DECISION: drop the matrix inset.
2. Ambient Hilbert space. Idea: one spin space over every configuration, with the norm as an integral.
   Drawing: a slice of C as a baseline with three configurations standing on it, a fibre line over
   each carrying the value Psi(X) as a two-component bar (labelled H_spin, never a plane), the
   component functions as one quiet plot beneath. 2D.
3. Physical Hilbert space. Idea: relabelling in the argument matches permuting the spin factors up to
   the statistics sign; parity is separate. Drawing: three water glyphs X, sigma X, E* X in a row with
   their spin boxes under the nuclei so the factors visibly follow the labels; the triad-on-sphere glyph
   (after Harter 5.1.1) for E* as inversion. 2D plus the orientation glyph.
4. Molecules, fragments, complexes. Idea: a grouping is a choice of boundaries on one configuration;
   inclusion is not accessibility. Drawing: methylamine as one fragment and as two (constructed 3D,
   gouache envelopes with heavy contour), KRb + KRb as a schematic channel scheme with K above and Rb
   below so no envelope crosses, channel groups as one nested-outline inset inside S.
5. Symmetry groups and versions. Idea: H merely rotates the reference; G adds resolved rearrangements;
   versions are G/H. Drawing: X0 with the mirror plane and the R_b axis constructed, the half-turn
   carrying X0 onto the ghost bX0 (verified: R_b = R_y(pi)); the bond interval [H, B] as a Hasse diagram
   with the chain as nested outlines (after 4.1.5) each holding the molecule and the extending
   operation; six versions as one solid reference and five ghosts with coset-leader labels at the sites.
6. Feasible configuration space. Idea: F is the thickened reference family carved out of C, with
   |S/G| disjoint images. Drawing: the short thick tube with a near-side cut-away wedge (band at
   mid-thickness, wall = Q, hatched cut faces, solid helical tu path, t path on the top edge, six gold
   versions, hidden ones on the bottom edge handled by leaders), a magnifying circle on the cut showing
   the normal slice; beneath, the version row as Newman projections down C-N along the t path; the
   twenty images as small solids inside the outline of C with the rest hatched as cut material.
   Disclosed: one radial thickness stands for 13 normal coordinates; SO(3) suppressed; the tube's
   topology is not a claim; assumptions are an embedded tubular neighbourhood and separation.
7. Symmetry species. Idea: three version states span A1 + E, one E subspace. Drawing: the coset
   triangle with t and b, and the three coefficient vectors as mode pictures (disc area = amplitude,
   fill = sign) after Harter 4.4.3; the E matrices as the only text. 2D.
8. Nuclear spin weights. Idea: each spatial species pairs with the spin species whose product contains
   the parity character; weights count copies. Drawing: the pairing matrix and the ruled table from
   the macros. DECISION: a 32-unit partitioned bar instead of the table.
9. Localized states. Idea: packet width sets the overlap, which sets the shares. Drawing: three gouache
   packets on the torsion circle with one contour each, the share curves to sigma/d = 1.2 with the
   displayed width marked, limits as open circles. 2D.
10. Position representation. Idea: near the family a configuration is r(X0(q) + sum Q e_a(q)) and a
    permutation-inversion acts on the chart by T_s. Drawing: the chart cut from a small tube glyph and
    unrolled with all six H-cells (uU split across the seam), q and q' with a straight arrow and the
    verified T_tu data; a magnifying circle for the normal slice with mass orthogonality; the
    reconstruction shown on one constructed molecule with declared arrow amplitude. DECISION: fold the
    slice into fig 6 and drop the reconstruction molecule.
11. Covering spaces and monodromy. Idea: twelve sheets in six version pairs over V; continuation
    multiplies on the right, so t-then-b and b-then-t end apart. Drawing: four oblique plates plus an
    ellipsis over the tinted base V, the two lifts drawn ending on different sheets with hidden parts
    dashed, the figure-eight loop on V (after Harter 5.4.2, disclosed as schematic), the two
    continuation rows beneath; X and bX as a ghost pair in a magnifying circle.
12. Momentum representation. Idea: the same states in modes; the seam shifts the torsional label to
    kappa = m + rho K; periodic torsion sorts into A1, A2 and E pairs. Drawing: the K rows with the
    K-ladder cone beside them (after 5.4.4, no energies), the species strip as drawn waveforms with
    level-rule typography (after 4.2.3), calligraphic J for the density. DECISION: drop the formula
    panel.
Guide: one nested carving drawing across the four parts, outside in, with one small emblem per part
and the exact part headings and section ranges; caption* so the counter is untouched. DECISION:
strip form without frames.

## 4. Repository layout, built to port to Overleaf
    figures/src/figNN-name.tex        standalone sources (pdflatex), one per plate plus fig00-guide
    figures/shared/figurestyle.sty    identity: colours, line styles, primitives
    figures/shared/primitives.tex     pics: atom, rod molecule, ghost molecule, thick tube, plates,
                                      triad sphere, zoom bubble, level rules, nested outlines
    figures/data/*.tex, *.dat         generated by compute/, committed (Overleaf needs no Python)
    figures/pdf/figNN-name.pdf        built plates, committed, for \includegraphics on Overleaf
    figures/preamble.tex              the lines the manuscript needs (\usepackage{standalone} or
                                      \graphicspath) and the \input or \includegraphics stanza per plate
    compute/, checks/, build.py       offline generation and verification (not needed on Overleaf)
    scaffold/outline.tex, proofsheet.tex
    docs/caption-notes.md             factual disclosures per plate, for the author's captions
    docs/verification.md              what was verified, how, and the residuals
    docs/author-decisions.md          the open decisions and the defaults taken
    docs/style.md                     the identity in one page
    README.md                         build, dependencies with tested versions, Overleaf porting
Overleaf porting: copy `figures/` and `docs/`; the manuscript includes each plate either by
`\input{figures/src/figNN}` with the standalone package (compiles on Overleaf, slower) or by
`\includegraphics{figures/pdf/figNN}` (instant). The integration patch is a separate, mechanical diff
against the author's manuscript.tex that adds only these lines. No shell escape, no network.

## 5. What I will do autonomously, in order, with gates
0. Fresh tree: new orphan branch `figures-v2` (or a new repository if the author creates one), with
   compute/, checks/, data/ carried over, everything else empty.
1. Primitives first: write primitives.tex and prove each pic on one primitives proof sheet at 150 mm,
   in colour, grayscale and CVD simulation. Gate: the sheet is sent to the author.
2. Prototypes 6, 5, 10 in the new idiom; render at scale; automated collision check (word boxes from
   `pdftotext -bbox`, flagged when any two overlap or any box leaves the plate); self-review; send the
   three renders. Gate: author comment; if none arrives, continue with the recorded defaults.
3. Remaining plates 1 to 4, 7 to 9, 11, 12; then the guide last. Same checks on each.
4. Clean-state build of every plate and the scaffold; proof sheet; numbering check (1 to 12, guide
   unnumbered); grayscale and CVD previews; collision check; a final look at every plate at 100 percent.
5. Documents: caption notes, verification record, author decisions, style page, README with tested
   versions; the integration patch; the source zip; commit and push.
Rules throughout: no prose or captions authored; no number typed by hand; one coordinator, at most two
bounded subagents, used only for isolated verification or a single plate; every delegated result is
rendered and inspected by the coordinator before it enters the tree.

## 6. Decisions the author owns (defaults in brackets, taken if unanswered)
- The five contract decisions above [keep the contracts as approved].
- Fresh branch in this repository or a new repository [branch `figures-v2`].
- Overleaf inclusion: inline standalone or built PDFs [both provided; PDFs recommended].
- Which Harter plates to cite beyond the adapted compositions [only the adapted ones].
- Ghost copies for the five non-reference versions, or six solid drawings [ghosts].
