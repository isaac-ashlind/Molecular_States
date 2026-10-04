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

Common to all plates: colours are the gouache set (ink, oxygen red, nitrogen blue, potassium green, rubidium
violet, gold for marked points and retained regions); blue solid lines are the path of t, red dashed lines the path
of tu (or b where stated), gold dots the version points; hatching marks cut or removed material only; molecules are
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
Coset triangle of G6/H with t and b. Mode pictures (disc area = amplitude, ink positive, white negative, a dash
for zero) are after Harter (1993), Fig. 4.4.3. Vectors v1, v2, v3 and the matrices of t and b on the E pair are
exact (verify_symbolic.py); C[G6/H] = A1 + E, no A2 (GAP).

## Figure 8, nuclear spin weights
Correlation diagram in the manner of Harter (1993), Fig. 4.2.3: spatial species as level rules, spin species as
bars of one unit per dimension (12, 4, 2 x 8 = 32), solid connections for even parity and dashed for odd. Copy
counts from the character computation (GAP, verify.py). chi_stat = +1 on G6.

## Figure 9, localized states
Three packets on the torsion circle at tau = 0, 2pi/3, 4pi/3: disc = one sigma, ring = 2 sigma, drawn with
sigma/d = 1/3 exactly; the share curves w_A1 = (1 + 2c)/3, w_E = 2(1 - c)/3 with c = exp(-d^2 / 8 sigma^2)
(figures/data/gaussian-shares.dat, exact in verify_symbolic.py); the drawn ratio is the marked point; limits as
open circles. Planar Gaussian illustration as in the manuscript.

## Figure 10, position representation
The band of figure 6 with two normal slices standing across it, cut at the seam and unrolled into the chart
[-pi, pi] x [-eta0, eta0] with the six H-cells (computed, uU split across the seam), the same six points and two
paths; T_tu carries q = (0, eta0) to q' = (-pi/3, -eta0) with A_tu = R_z(pi) and M_tu = diag(1, -1) (verified
identically); the normal slice at q is the cut face of figure 6; e_-(q) drawn on one computed molecule with
arrows 3 x the displacement Q_- = 0.42 A (`\arrowGain`, `\dispAmp`). Mass orthogonality holds in the mass metric,
not on the page.

## Figure 11, covering spaces and monodromy
Twelve sheets over V in six version pairs; lifts of the two loops drawn climbing to different sheets: a closed
path downstairs whose lift does not close, in the spirit of Harter (1993), Fig. 5.4.2. The relations t^3 = b^2 = E,
btb = t^-1, bt = t^2 b, the Schreier graph, the residual 0.76 A (no rotation carries X to bX) and the free-action
margin 2.30 A are computed. The plates are schematic; V is a free domain by assumption.

## Figure 12, momentum representation
K ladder as latitude rings on the sphere |J| = sqrt(J(J+1)) for J = 2, after Harter (1993), Fig. 5.5.3; no
energies. The (m, K) lattice in the periodic frame and sheared in the twisted frame with kappa = m + rho K,
rho = 0.3 illustrative (declared). Periodic torsion modes as radial plots on the torsion circle of figure 9;
species from U_t f_m = e^{-2 pi i m/3} f_m, U_b f_m = f_{-m} (computed); level rules after Harter (1993),
Fig. 4.2.3.
