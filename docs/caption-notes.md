# Caption notes (disclosures for the author's captions)

No final captions are written here. For each plate: what is drawn, what is computed (and by which check), what
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
is a surface and at full strength for rubidium; gold for retained regions and for potassium; white cladding and
hydrogen); actions are ink, solid for t, dashed for u or tu, dotted for a starred operation (b, E*); red dots are
the version points; hatching marks cut or removed material only; molecules are
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
15 degrees so that the projection line from H to [X] crosses no other version point.

## Figure 1, configuration space
Water atlas: 25 glyphs computed from (delta, theta) with r1 = l + delta, r2 = l - delta, delta in {-l/2, -l/4, 0,
l/4, l/2}, theta in {60, 90, 120, 150, 180} degrees, each centred on its own mass centre (no mark drawn, 0.06 A from
the oxygen in the upper-left shape, `\waterOxygenOffset`); dotted loci delta = 0 and theta = pi, their intersection
circled. Commuting square: R = R_z(40 deg) on the left, P_sigma with sigma = (12) on the right; the matrix X and
X m = 0 from the data (`\waterX..`, `\massH`, `\massO`). Verified in checks/verify.py (centring, P_{st} = P_t P_s).

## Figure 2, ambient Hilbert space
A slice of C (delta = 0, theta from 80 to 180 degrees) as a baseline with a comb of fibres, one spin space over every
configuration; on each fibre the value Psi(X) as its two components a (teal) and b (gold) in the basis
|up down> omega, |down up> omega; three configurations X1, X2, X3 standing on the slice; the component functions
plotted beneath on the same theta scale with the norm density dotted and the three sampled values as dots. The
state drawn is physical on this slice by the author's choice: rotation-invariant and exchange-antisymmetric, so
b = -a (the proton singlet times f(theta), f an illustrative Gaussian, declared). Nothing in the ambient space
forces this; the caption should say the plate shows a physical representative, not a generic one.

## Figure 3, physical Hilbert space
Top row: X with spin arrows beside the nuclei (up teal on 1, down gold on 2, the oxygen spinless); sigma X = X P_sigma
with the digits moved and the spins left with the places; -X P_sigma = E* sigma X with every column sent through the
origin. The two operations as arrows carrying their sign: sigma = (12) with chi_stat(sigma) = -1 (solid), E* with
chi_pm(E*) = +-1 (dotted, a starred operation). Bottom row: the spin value itself, two coloured arrows in a bracket
in slot order (1, 2): Psi(X); Psi(X P_sigma) = chi_stat(sigma) P_sigma Psi(X), the arrows permuted and the sign in
front; Psi(-X P_sigma) = chi_pm Psi(X P_sigma), the parity sign in front. Water is planar, so -X P_sigma is also a
rotated copy of X (caption, not plate). Configurations computed (compute/geometry.py).

## Figure 4, molecules, fragments, complexes
The reaction plate. One set of four nuclei (potassium 1, 2; rubidium 3, 4) grouped three ways along the reaction,
KRb + KRb (G_in = <(12)(34), E*>), the complex K2Rb2 (S = <(12), (34), E*>), K2 + Rb2 (G_K = <(12), E*>,
G_Rb = <(34), E*>); a line joins grouped nuclei; positions schematic (declared). The two product channels are two
strands, gold for potassium and teal for rubidium, side by side from the reactants into the complex and split at
the products; the same strands on the lattice of the sixteen subgroups of S (congruent to C_2^3, computed in
checks/verify.g), up from G_in to S and down to G_K and G_Rb; the common core <E*> below, E at the bottom; the
lattice shows inclusion, not accessibility. Named nodes gold (G_K), teal (G_Rb), half and half (G_in).

## Figure 5, symmetry groups and versions
Row A: methylamine as two fragments, methyl and amino (dashed hulls); X0 with its half-turn axis and the version
equation b = (23)(45)* = R_y(pi), bX0 = X0 P_b = R_b X0, so H = {E, b}. Row B: the 36 subgroups between H and S
(grey) carved to the 10 that keep the bonds (black), the chain H < G6 < G12 < B in gold; at the right the bond
interval [H, B] as a Hasse diagram, every node written by its generators, each cover drawn in the line grammar of
the generator it adjoins (solid t, dashed u, dotted a starred one: E*, u*, t*), the same chain as a gold halo.
Row C: the six versions G12/H, the same positions with the digits moved, position j carrying label g(j). All
lattice data, generator forms and cover labels computed (compute/, checks/verify.g).

## Figure 6, feasible configuration space
The solid F as a short thick tube cut open: band at mid-thickness, tau around, eta along the height, Q across the
wall; one radial thickness stands for the 13 normal coordinates; SO(3) suppressed; the tube's topology is not a
claim about F; assumed: an embedded tubular neighbourhood and separation of the images. Version positions, the
t and tu paths and the eclipsed point are the computed family positions. Newman projections are computed from the
family (digits are column labels). Twenty images since |S/G12| = 20 (GAP).

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
bundle it is joined to: 12, 4, 8 (even A1, A2, E) and 4, 12, 8 (odd). chi_stat = +1 on G6. Counts from the character
computation (GAP, verify.py); the same counts arise as the product of the methyl multiplets (4 A1 + 2 E) and the
amino multiplets (3 + 1) (checks/verify_spin.py), which the caption may mention.

## Figure 9, localized states
The torsion circle unrolled to [-pi, pi]; three periodic Gaussians at tau = 0, 2pi/3, 4pi/3 in three rows,
sigma/d = 1/6, 1/3 (the ratio drawn elsewhere), 2/3; the overlap c = <g0, g1> = exp(-d^2 / 8 sigma^2) is the darker
gold; the bar at the right of each row is the A1 share w_A1 = (1 + 2c)/3 (filled) and the E share w_E = 2(1 - c)/3
(open); the small chart marks the three rows on the share curves (figures/data/gaussian-shares.dat, exact in
verify_symbolic.py). Planar Gaussian illustration as in the manuscript.

## Figure 10, position representation
The band of figure 6 (the reference family X0(tau, eta) at mid-thickness) with its six version points and the two
paths, cut at the seam tau = +-pi (cut marks) and unrolled into the chart [-pi, pi] x [-eta0, eta0]; the six
H-cells computed (chartCells), the reference cell U in pale gold, the same six red points and the two paths (t
solid, tu dashed); eta = 0 is the planar amino locus. The zoom at q = H is the normal slice, the hatched cut face
of figure 6, with the band's eta-line (teal) and the coordinate Q across the wall. Below, the two normal
directions at q, e_-(q) (the amino twist) and e_+(q) (the amino rock), drawn on the molecule with arrows 3x the
displacement (mol3d-mla-ref-eminus, -eplus; mass-orthonormal, orthogonal to the family and to the rotations in
the mass metric). The formulas X0(q) P_tu = A_tu X0(q'), T_tu and M_tu = diag(1, -1) are in the text, not the plate.

## Figure 11, covering spaces and monodromy
Over a free domain V, the version pairs {gX, gbX} drawn side by side: the X-side column at the left headed by X,
the bX-side at the right headed by bX (two configurations over the one point [X]; no rotation carries one to the
other, best-fit residual 0.76 A; the action is free near X0, nearest image 2.30 A). The loop l_t (solid) lifts by
climbing a column, l_b (dotted) by hopping across; from X, t then b ends at tbX one level up on the right, b then t
ends at btX = t^2 bX two levels up on the right (gold rings): the lifts of the same loops end apart. At the right the
six t,b-sheets as a graph: t turns the outer triangle one way and the inner the other (btb = t^-1). The u-half of
the twelve sheets is alike. Loop figure after Harter (1993), Fig. 5.4.2 in spirit.

## Figure 12, momentum representation
K ladder as latitude rings on the sphere |J| = sqrt(J(J+1)) for J = 2, after Harter (1993), Fig. 5.5.3; no
energies. The (m, K) lattice in the periodic frame and sheared in the twisted frame with kappa = m + rho K,
rho = 0.3 illustrative (declared). Periodic torsion modes as radial plots on the torsion circle of figure 9;
species from U_t f_m = e^{-2 pi i m/3} f_m, U_b f_m = f_{-m} (computed); level rules after Harter (1993),
Fig. 4.2.3.
