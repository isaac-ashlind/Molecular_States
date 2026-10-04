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

## Figure 1, configuration space
Water atlas: 25 glyphs computed from (delta, theta) with r1 = l + delta, r2 = l - delta, delta in {-l/2, -l/4, 0,
l/4, l/2}, theta in {60, 90, 120, 150, 180} degrees, each centred on its own mass centre (gold cross, 0.06 A from
the oxygen in the upper-left shape, `\waterOxygenOffset`); dotted loci delta = 0 and theta = pi, their intersection
circled. Commuting square: R = R_z(40 deg) on the left, P_sigma with sigma = (12) on the right; the matrix X and
X m = 0 from the data (`\waterX..`, `\massH`, `\massO`). Verified in checks/verify.py (centring, P_{st} = P_t P_s).

## Figure 2, ambient Hilbert space
A slice of C (delta = 0, theta from 80 to 180 degrees) as a baseline with a comb of fibres; on each fibre the two
components a(X), b(X) of one continuous representative (illustrative functions, declared); three configurations at
theta = 90, 110, 140 degrees stand on the slice; the plot beneath uses the same theta scale and the same functions
(figures/data/component-functions.dat). The norm is the integral of |a|^2 + |b|^2. Nothing here imposes exchange.

## Figure 3, physical Hilbert space
X, sigma X = X P_sigma and E* X = -X as computed configurations; the spin slots in label order with leaders to the
nuclei carrying each digit; the statistics sign and the parity rule as the only text. Water is planar, so -X is
also a rotated copy (stated on the plate); the inversion triad on a sphere is after Harter (1993), Fig. 5.1.1.

## Figure 4, molecules, fragments, complexes
Methylamine once with one fragment boundary and once with two (fine dashed contours around the computed
positions). K2Rb2: positions schematic (declared), potassium above and rubidium below, a thin line joins grouped
nuclei; the three channel groups, their orders, their common core <E*>, the fact that any two generate S, and the
sixteen subgroups of S with their covers are computed (compute/groups.py, checks/verify.g).

## Figure 5, symmetry groups and versions
X0 with the mirror plane through H1, C, N and the half-turn axis; bX0 = X0 P_b drawn solid with digits, equal to
R_y(pi) X0 to 1e-15 (verify.py, verify_symbolic.py). The 36 subgroups between H and S (grey) with the 10 of the
bond interval [H, B] (black) and the chain H < G6 < G12 < B (blue), the same chain in the Hasse diagram of [H, B]
(10 subgroups, 17 covers); all computed and cross-checked in GAP. Six versions as the same positions with the
digits moved by g (position j carries label g(j)).

## Figure 6, feasible configuration space
The solid F as a short thick tube cut open: band at mid-thickness, tau around, eta along the height, Q across the
wall; one radial thickness stands for the 13 normal coordinates; SO(3) suppressed; the tube's topology is not a
claim about F; assumed: an embedded tubular neighbourhood and separation of the images. Version positions, the
t and tu paths and the eclipsed point are the computed family positions. Newman projections are computed from the
family (digits are column labels). Twenty images since |S/G12| = 20 (GAP).

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
