# Caption notes (disclosures for the author's captions)

Symbols below are written in plain ASCII. In the plates and the manuscript they follow the approved notation
(2026-10-04): configuration matrices X, X0 and the action matrices P_s, R, A_s, M_s are bold; the spatial operators
U_s and the projector P_Gamma are sans serif; the umbrella coordinate is iota, the shape tuple a, the normal
coordinates q, the packet width Delta.

The manuscript carries minimal captions (one to three sentences each, 2026-10-04, author's request); these notes
remain the disclosures behind them. For each plate: what is drawn, what is computed (and by which check), what
is declared schematic or illustrative, and which compositions adapt a plate of W. G. Harter, *Principles of
Symmetry, Dynamics, and Spectroscopy* (Wiley, 1993), to be cited as "after Harter (1993), Fig. x.y.z".

```bibtex
@book{Harter1993,
  author    = {Harter, William G.},
  title     = {Principles of Symmetry, Dynamics, and Spectroscopy},
  publisher = {Wiley},
  address   = {New York},
  year      = {1993}
}
```

Common to all plates: colours are the anchored gouache set (ink; red for oxygen on the water plates and for the
version points on the methylamine plates; blue for nitrogen; teal as a pale tint for the reference family where it
is a surface and at full strength for rubidium; gold for retained regions and for potassium; white vibrations and
hydrogen); actions are ink, solid for t, dashed for u or tu, dotted for a starred operation (b, E*); red dots are
the version points; hatching marks cut or removed material only; every line is ink, varied by weight and style,
and grey is kept for greyed-out objects (the cover's lattice, the inactive parts of the lattices of figures 4 and 5,
the undisplaced shape under the displaced one in figure 10, the hidden version points); molecules are
drawn from the computed configurations at 0.95 cm per angstrom in one oblique projection (camera azimuth 50 deg,
elevation 25 deg); digits inside hydrogen discs are the column labels 1 to 5.

## Guide (unnumbered)
One scene: the configuration space as a hatched slab with the hole the carving leaves, the feasible solid lifted
out of it with its band, the two paths and the six version points (hidden ones faded), the covering sheets lifted
above one version, one configuration as a point of the space. Schematic throughout; the solid is the figure 6
object. Set with `\caption*` so the figure counter is untouched.

The six version points on the core sheet are the vertices of a trigonal antiprism; their symmetry group within the
cylinder's symmetries is of order 12 and isomorphic to G12 (checks/verify_antiprism.py), so the decorated solid
carries exactly the molecular symmetry group, which the caption may say. The solid's decoration is drawn turned
15 degrees so that the projection line from H to [X] crosses no other version point. The four headings carry the
roman numerals of the four sections (the former parts) with their subsection ranges. At the top left, under heading
III, the bond interval [H, B] of figure 5 as a small Hasse diagram in grey (ink thinned), covers in the line grammar
(solid t, dashed u, dotted a starred generator), the vertices H, B and G (= G_12, the regime) named; no colour. The
base patch of the sheets is pale gold, a chart neighbourhood in the quotient (gold = a retained region, here the
neighbourhood carrying a chart), with the class of the reference as a red point; it is unnamed (the quotient is named
on the projection, F -> F/G, the same G as the lattice vertex), and the dotted projection runs unbroken from the red
point on the sheet to the point on the patch. The molecule is drawn mass-centred, and the dotted drop from it to its
point of C is the vertical through its centre of mass (the same dotted line as the projection from H to the patch).

## Figure 1, configuration space
Water atlas: 25 glyphs computed from (delta, theta) with r1 = l + delta, r2 = l - delta, delta in {-l/2, -l/4, 0,
l/4, l/2}, theta in {60, 90, 120, 150, 180} degrees, each centred on its own mass centre (no mark drawn, 0.06 A from
the oxygen in the upper-left shape, `\waterOxygenOffset`); the loci delta = 0 (equal bonds) and theta = pi (linear) are
the atlas's own coordinate lines, drawn as thin solid ink (thinner than the axes; dashes and dots are kept for actions),
not a schematic boundary. Commuting square: the rotation matrix R = R_z(40 deg) on the left, the column matrix P_sigma with sigma = (12) on the right; the
matrix X and X m = 0 from the data (`\waterX..`, `\massH`, `\massO`). Verified in checks/verify.py (centring, P_{st} = P_t P_s).

## Figure 2, ambient Hilbert space
A slice of C (delta = 0, theta from 80 to 180 degrees) as a baseline with a comb of fibres, one spin space over every
configuration; on each fibre the value Psi(X) as its two components a (teal) and b (gold) in the basis
|up down> omega, |down up> omega (the text's xi and zeta; zeta replaces the earlier eta, which is the amino
coordinate elsewhere); three configurations X1, X2, X3 standing on the slice, a sample of the one-parameter family along the slice, not a
special set; the component functions plotted beneath on the same theta scale with the norm density |a|^2 + |b|^2 as
a thin ink curve (no marks: each sample stands directly above its density value because the scales agree). The
state drawn is physical on this slice by the author's choice: rotation-invariant and exchange-antisymmetric, so
b = -a (the proton singlet times f(theta), f an illustrative Gaussian, declared). Nothing in the ambient space
forces this; the caption should say the plate shows a physical representative, not a generic one.

## Figure 3, physical Hilbert space
Top row: X with spin arrows beside the nuclei (up teal on 1, down gold on 2, the oxygen spinless); sigma X = X P_sigma
with the digits moved and the spins left with the places; -X P_sigma = E* sigma X with every column sent through the
origin. The two operations as arrows carrying their sign: sigma = (12) with chi_stat(sigma) = -1 (solid), E* with
chi_pm(E*) = +-1 (dotted, a starred operation). Bottom row: the spin value itself, two coloured arrows in a bracket
in slot order (1, 2): Psi(X); Psi(X P_sigma) = chi_stat(sigma) sigma.Psi(X), the arrows permuted (the spin action of
sigma, as the text writes it) and the sign in front; Psi(-X P_sigma) = chi_pm Psi(X P_sigma), the parity sign in front. The
plate fixes one parity sector at a time (the sign chi_pm is a choice of sector, written +- on the plate). Water is
planar, so -X P_sigma is also a rotated copy of X (caption, not plate). Configurations computed (compute/geometry.py).

## Figure 4, molecules, fragments, complexes
The reaction plate. One set of four nuclei (potassium 1, 2; rubidium 3, 4) grouped three ways along the reaction,
KRb + KRb (G_in = <(12)(34), E*>), the complex K2Rb2 (S = <(12), (34), E*>), K2 + Rb2 (G_K = <(12), E*>,
G_Rb = <(34), E*>); a line joins grouped nuclei; positions schematic (declared). The two product channels are two
strands, gold for potassium and teal for rubidium, side by side from the reactants into the complex and split at
the products; the same strands on the lattice of the sixteen subgroups of S (congruent to C_2^3, computed in
checks/verify.g; the lattice is complete, every subgroup is a node, the unnamed ones small and grey), up from G_in to S
and down to G_K and G_Rb; the common core <E*> below, E at the bottom; the lattice shows inclusion, not accessibility. Named nodes gold (G_K), teal (G_Rb), half and half (G_in).

## Figure 5, symmetry groups and versions
Row A: methylamine as two fragments, methyl and amino (dashed hulls); X0 with its half-turn axis and the version
equation b = (23)(45)* = R_y(pi), bX0 = X0 P_b = R_b X0, so H = {E, b}. Row B: the 36 subgroups between H and S
(grey) carved to the 10 that keep the bonds (black), the chain H < G6 < G12 < B in gold; at the right the bond
interval [H, B] as a Hasse diagram, every node written by its generators, each cover drawn in the line grammar of
the generator it adjoins (solid t, dashed u, dotted a starred one: E*, u*, t*), the same chain as a gold halo.
Row C: the six versions G12/H, the same positions with the digits moved, position j carrying label g(j); each glyph is
one representative of its rotational-orbit version (the orbit under H = {E, b} is the rotation by R_b). All
lattice data, generator forms and cover labels computed (compute/, checks/verify.g).

## Figure 6, feasible configuration space
The vibrations are the displacements mass-orthogonal to the rotations and to the family's tangent plane: at each
shape the tangent space of C splits into 3 rotational directions, 2 along the family (tau, iota) and the 13 normal
ones; the wall's thickness is the reach of those 13, and the zoom's q arrow runs along that normal direction from
the sheet's edge across the wall.
The solid F as a short thick tube cut open: band at mid-thickness, tau around, iota along the height, q across the
wall; one radial thickness stands for the 13 normal coordinates; SO(3) suppressed; the tube's topology is not a
claim about F; assumed: an embedded tubular neighbourhood and separation of the images. Version positions, the
t and tu paths and the eclipsed point are the computed family positions. Newman projections are computed from the
family (digits are column labels). Twenty images since |S/G12| = 20 (GAP); the active image is a plain white miniature
of the solid, the others grey.

The Newman row shows the two steps between versions as one arrow each from H: tu (dashed) to tuH and t (solid) to
tH. Each is a single relabelling that sends one version straight to another. The eclipsed shape (tau = pi/3) stands
apart at the right as a reference: it is the midpoint in tau of the continuous path that realizes t on the sheet,
not a version, and t does not pass through it in steps. On the solid the solid path from H to tH is that continuous
path, with the eclipsed point marked white on it.

In the eclipsed Newman projection (tau = pi/3) the back set is drawn turned 16 degrees on the page so that both
sets of hydrogens show, the usual drawing convention (declared; the geometry itself is exactly eclipsed).

## Figure 7, symmetry species
The three version sites on the torsion circle (tau = 0, 2pi/3, 4pi/3); t turns them by a third, b is the mirror
through H (dotted line, the starred generator). Mode pictures (disc area = amplitude, ink positive, white negative,
a dash for zero) are after Harter (1993), Fig. 4.4.3 in spirit; the grid shows each mode v, its image t v and its
image b v; the plane beneath shows the three images of v2 at a third turn with the mirror axis. Vectors v1 =
(1,1,1)/sqrt3, v2 = (2,-1,-1)/sqrt6, v3 = (0,1,-1)/sqrt2 and the actions t(a,b,c) = (c,a,b), b(a,b,c) = (a,c,b)
are exact (verify_symbolic.py); C[G6/H] = A1 + E, no A2 (GAP).

## Figure 8, nuclear spin weights
Correlation diagram in the manner of Harter (1993), Fig. 4.2.3: the proton spin space H_spin = (C^2)^(x5) in the
centre as bundles of hairlines, one line per state, 12 A1, 4 A2 and 8 E pairs (32 in all); even-parity spatial
species (chi_+ = A1) at the left joined straight across to the spin species of the same name; odd-parity species
(chi_- = A2) at the right joined with A1 and A2 swapped and E kept; a level's weight is the number of lines in the
bundle it is joined to: 12, 4, 8 (even A1, A2, E) and 4, 12, 8 (odd); these are multiplicities (copies of a species),
not dimensions, so the eight E pairs fill sixteen of the thirty-two dimensions. chi_stat = +1 on G6 (t is a 3-cycle
and the permutation part of b is a double transposition). Counts from the character
computation (GAP, verify.py); the same counts arise as the product of the methyl multiplets (4 A1 + 2 E) and the
amino multiplets (3 + 1) (checks/verify_spin.py), which the caption may mention.

## Figure 9, localized states
The torsion circle unrolled to [-pi, pi]; three periodic Gaussians at tau = 0, 2pi/3, 4pi/3 in three rows
of Delta/d = 1/6, 1/3 (the ratio drawn elsewhere), 2/3 (Delta is the packet width, the density standard deviation per
coordinate; sigma stays the permutation symbol); the overlap c = <g0, g1> = exp(-d^2 / 8 Delta^2) is the
gold; the bar at the right of each row is the A1 share w_A1 = (1 + 2c)/3 (filled) and the E share w_E = 2(1 - c)/3
(open); the small chart marks the three rows on the share curves, the A1 share heavy and the E share thin, both ink,
with the three ratios as ticks of the Delta/d axis (figures/data/gaussian-shares.dat, exact in verify_symbolic.py).
The packets are illustrative torsional profiles drawn on the unrolled circle; the overlap and the shares use the
planar Gaussian formula of the manuscript with the packet width Delta, not a torsional potential, and the caption
should say so.

## Figure 10, position representation
Top row, left to right: the solid F of figure 6 (cut open, window -62 to 8 degrees, the decoration unturned so that the
seam sits at the back), stripped of its vibrations to the core sheet (a thin teal tube of the sheet's own thickness, the
reference family X0(a), an open cylinder in a = (tau, iota), cut at the seam tau = +-pi, cut marks; the small solid
carries the version points but not the paths, which are read on the sheet and the chart), unrolled into the chart [-pi, pi] x
[-iota0, iota0] with the six H-cells computed (chartCells), the reference cell U in pale gold, the same six red
version points and the two paths (t solid, tu dashed); iota = 0 is the planar amino locus. Middle row: over a point
a of the chart sit every normal displacement q and every orientation r, drawn as the reconstruction: the shape X0(a),
the shape displaced along the normal frame X0(a) + sum q e(a) (amplitude 0.42 A along e_-, illustrative, declared;
its undisplaced shape beneath as faded solid outlines), and the result rotated, phi(r,a,q) = R(r)(...) (rotation about (0.3, 1, 0.25) by
55 degrees, illustrative, declared). Bottom row: the two normal directions at a, e_-(a) (the amino twist, hydrogens 4
and 5 against each other along the axis) and e_+(a) (the C-N stretch, carbon against nitrogen along the axis), on
the molecule seen from the side (camera azimuth 75, elevation 18, so the C-N axis lies in the page and no methyl
hydrogen sits on the line of sight through the carbon); each direction's arrows are scaled so that its largest arrow
is 1.1 A on the page (the direction is the content; the mass-orthonormal stretch moves the heavy atoms far less than
the twist moves the hydrogens), declared (mol3d-mla-ref-eminus-side, -eplus-side; mass-orthonormal, orthogonal to
the family and to the rotations in the mass metric). The formulas X0(a) P_tu = A_tu X0(a'), T_tu and
M_tu = diag(1, -1) are in the text, not the plate.

## Figure 11, covering spaces and monodromy
Over a patch of the quotient F/G, drawn the size of one sheet and pale gold like the chart cell of figure 10 (each
sheet maps onto it), the version pairs {gX, gbX}
side by side: the X-side column at the left headed by X, the bX-side at the right headed by bX (two configurations over
the one point [X], a red point like the version points; no rotation carries one to the
other, best-fit residual 0.76 A; the action is free near X0, nearest image 2.30 A). The loop l_t (solid) lifts by
climbing a column, l_b (dotted) by hopping across; from X, t then b ends at tbX one level up on the right, b then t
ends at btX = t^2 bX two levels up on the right: the lifts of the same loops end apart. The rule behind both: the group
acts on the left, so the lift of l_t that starts at gX ends at gtX, and from bX that is btX = t^2 bX because btb = t^-1;
"t then b" names the order of the two moves from X. The six sheets drawn are the G6 block of the twelve (one per
element of G6 over the point [X], rows named by the cosets of H); the u-half is alike. At the right the six sheets as
a graph: the solid arrows are the lifts of l_t, which turn the outer triangle one way and the inner the other; the
dotted edges are the lifts of l_b. The projections from [X] to X and bX are dotted like the cover's. Loop figure after
Harter (1993), Fig. 5.4.2 in spirit.

## Figure 12, momentum representation
K ladder as latitude rings on the sphere |J| = sqrt(J(J+1)) hbar for J = 2 (the plate writes the hbar), after Harter
(1993), Fig. 5.5.3; no energies, the rings are the values of K only. The (m, K) lattice in the periodic frame and sheared in the twisted frame with kappa = m + rho K and
rho = 1/2, the symmetric frame of Mellor, Yurchenko, Mant and Jensen (Symmetry 11, 862, 2019), so the twisted
labels are half-integers; the zero line and the gold line kappa = rho K lie under the lattice dots. Periodic torsion
modes as radial plots (heavier ink) on the torsion circle of figure 9 (thin ink, the three version points at the
suite's radius);
species from U_t f_m = e^{-2 pi i m/3} f_m, U_b f_m = f_{-m} (computed); level rules after Harter (1993),
Fig. 4.2.3.

## Rigid recovery (fig13, unnumbered, the three closing pages): methane through the pipeline, against Albert et al.
What the reader sees at a squint: a labelled tetrahedron with three axes; a circle holding a scatter of red dots, a
small wire octahedron and one arrow; a stack of three flat sheets over a base patch with a loop and a lift. The ball
panel has its own camera (make_data: CAM_BALL, azimuth 108, elevation 21, chosen by scan to sit 27 degrees, the
maximum, from every symmetry axis of the cube, so the twelve sites and the six cell vertices all separate and bond 1
faces the viewer); the molecule panel has the molecule's (CAM_MOL, azimuth 18, elevation 59, chosen by scan so that
every hydrogen disc clears the carbon and every other hydrogen by 0.61 A on the page with every C2 axis at least 26
degrees off the view). The two panels share the body frame, not the camera. Left: the reference X0 in its C2 frame
(geometry.methane_c2: hydrogens on the body diagonals, C-H 1.087 A) drawn with the suite's molecule primitive
(mol3d-ch4-c2: shaded atoms, rods, the digits inside the hydrogens), the axes x, y, z through the midpoints of
opposite edges as fine lines behind it.
Middle: the orientation space SO(3) as the axis-angle ball (direction the axis, distance the angle, radius pi drawn
2 cm; the skin is the half-turns, each point its own antipode): the twelve symmetry rotations of T as red dots, the
centre, eight third-turns (a tetrahedron through the bonds at +2pi/3, one through the faces at -2pi/3) joined by the
thin edges of the cube they span, and three half-turns on the skin at the crossings of the great circles, each drawn
once at its visible representative; the three great circles of the coordinate planes give the depth; one cell of
the twelve, the orientations nearer the centre than any other dot, as a wire octahedron (an octahedron in Rodrigues
coordinates, |x|+|y|+|z| <= 1; drawn with flat faces, the true faces bulge) with the six quarter-turns, the starred
4-cycles, as ink rings at its vertices; the lift of the
third-turn loop as a solid arrow from the centre through a face to the dot beyond (r -> r R_g). No molecule inside
the ball and no point drawn twice (author, 2026-10-05). Right: the covering in figure 11's grammar: a patch of
SO(3)/T (the cell with its faces glued) with the loop l_g, the sheets r, r R_g, r R_g^2 of the twelve over it, the
lift of l_g from the dot on sheet r to the dot on sheet r R_g, one dotted projection from the base point to its
sheet; the sheets are the cartoon of figure 11, not a computed surface. Line vocabulary (author, 2026-10-05): solid
for what faces the viewer, dashed for what lies behind; dotted only for the projection. Colours untangle the nested
surfaces (author): the sphere teal, the orientation space taking over the core sheet's role; the cell gold, the
chart patch as in figures 10 and 11, and the base patch at right in figure 11's gold tint; the cube in thin ink; the
dots red; the lift in ink. This is the honest form of Albert et al.'s
Fig. 3(c): their fibre over a position is these twelve sheets (the regular representation of T), and their
monodromy action is the lift moving from one sheet to the next. The red dots carry the version-dot role of figure 10,
handed to orientations. Every coordinate is emitted by make_data (numbers.tex, the ball*
macros) and the geometry is recomputed in checks/verify_methane.py: all 24 relabellings are proper rotations
closing as the cube group O, the even twelve as T; the identity's Voronoi cell in the bi-invariant metric is the
octahedron (3000 random rotations); its vertices are the quarter-turns; the face x+y+z=1 is glued to x+y+z=-1 by
R -> R R_g^{-1} with the vertex R_x(pi/2) landing at R_z(-pi/2), a third of a twist; the loop about bond 1 returns
X0 with 2, 3, 4 cycled and R_g cycles x, y, z. The three pages carry the pipeline (sections 1 to 10 instantiated on
rigid methane; the vibrational species A1 + E + 2 T2 from explicit 15 x 15 matrices; the packet shares
(1+8c3+3c2)/12, (1-4c3+3c2)/12, 3(1-c2)/4; the J ladder D^J restricted to H with the physical states by parity,
J = 0..6) and the projection onto Albert et al. (Tables A to D: the classes of Td(M); the species with their isomers,
their Table II T block; the J ladder; the dictionary with their equation and example numbers). Declared: no
energies; parity (the starred half of H) and the vibrations (the disc D^9 under M_h) are what their formalism sets
aside; the labels 1E, 2E are exchanged together by the other reading of the transport (footnote).
