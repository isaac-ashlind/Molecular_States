# Verification

Every number and structural claim drawn in the plates is computed in `compute/` and checked three ways.

| Layer | File | What it establishes |
|---|---|---|
| exact and floating arithmetic | `checks/verify.py` | group orders and cosets; [H,B] has 10 subgroups and 17 covers; [H,S] has 36 subgroups and 73 covers; relations btb = t^-1, bt = t^2 b; spin weights and C[G6/H] multiplicities; P_{st} = P_t P_s; the E-pair matrices; torsion species; Gaussian shares and limits; centring to 1e-12; b = R_y(pi) to 1e-15; family actions to 1e-15; min |gX0 - X0| = 2.30 A; emitted data files agree with recomputation; KRb groups |
| symbolic | `checks/verify_symbolic.py` (SymPy) | the four family actions hold identically in (tau, eta); X m = 0 identically; c = exp(-d^2/8 Delta^2) and the shares exactly; the E matrices exactly; the seam relation kappa = m + rho K exactly |
| independent algebra | `checks/verify.g` (GAP) | orders 2, 6, 12, 24, 240; cosets 3, 6, 20, 120; [H,B] 10/17 with orders [2,4,4,4,6,8,12,12,12,24]; three order-12 extensions, any two generate B; [H,S] 36/73; class sizes; spin weights 12/4/8 (even) and 4/12/8 (odd); C[G6/H] = A1 + E; KRb: |S| = 8, 16 subgroups, channel orders 4 |
| spin structure | `checks/verify_spin.py` (numpy) | the five proton spins as methyl times amino: 4 A1 + 2 E, 3 + 1, and 12 A1 + 4 A2 + 8 E under G6 (below) |
| band symmetry | `checks/verify_antiprism.py` | the decorated band of the cover and figure 6 has exactly G12 as its symmetry (below) |
| rigid methane | `checks/verify_methane.py` (numpy) | every number and construction on the closing pages (below) |
| text | `checks/collisions.py` | overlapping word boxes per plate (pdftotext -bbox); remaining reports are glyph pairs of one symbol (the ≅ sign, roots in fractions) and the eclipsed digits of figure 6 |

`python3 build.py` runs all of them before compiling; a failing check stops the build. Declared schematic or illustrative: K2Rb2 positions (figure 4), the component functions a, b
(figure 2), rho = 0.3 (figure 12), the tube as a picture of F (figure 6), the sheets (figure 11).

## Spin-weight structure (figure 8, added 2026-10-04)

`checks/verify_spin.py` builds the permutation matrices on the 2^n spin basis and decomposes by characters:
the three methyl spins under S3 give 4 A1 + 2 E (the quartet and the two doublets); the two amino spins under the
swap give 3 symmetric + 1 antisymmetric (triplet, singlet); the five together under G6 = <(123), (23)(45)> give
12 A1 + 4 A2 + 8 E, equal to the product (4 A1 + 2 E) x (3 A1 + 1 A2). This is the mosaic of figure 8: 8 methyl
states across by 4 amino states down, four blocks of 12, 4, 12 (6 E pairs) and 4 (2 E pairs) cells. It agrees
with the GAP character computation (checks/verify.g) that produced the counts in figures/data/numbers.tex.

## The decorated band as a trigonal antiprism (cover, figure 6; added 2026-10-04)

`checks/verify_antiprism.py`: G12 = <b, t, u> modelled as signed permutations has order 12 with element orders
(1, 7, 2, 2 of orders 1, 2, 3, 6), the profile of D6 = S3 x C2; it acts faithfully on the six cosets G12/H, H = <b>,
i.e. on the six version points; and the stabiliser of the six points {(0,+), (2pi/3,+), (4pi/3,+), (pi,-), (5pi/3,-),
(pi/3,-)} inside the symmetries of the (tau, eta) cylinder (rotations by multiples of pi/3, tau -> -tau, eta -> -eta;
order 24) has order 12 with the same profile. So the band decorated by the orbit of the reference has exactly the
molecular symmetry group as its symmetry: t is the third turn, u the inversion through the centre (tau + pi, -eta),
tu the sixfold rotation-reflection, b a vertical mirror; the six points are the vertices of a trigonal antiprism
(point group D3d, isomorphic to G12). The undecorated band has the continuous symmetry of a cylinder; the orbit cuts
it down to G12 and to nothing less, because the reference X0 is generic (its stabiliser in G12 is only H).

## Rigid methane (the closing pages; added 2026-10-05)

`checks/verify_methane.py` recomputes, from explicit permutation, rotation and displacement matrices, every number on
the closing pages and compares it with what `compute/make_data.py` emits (`figures/data/summary.json`, keys
`methane_rigid`, `methane_loop`):
- every relabelling of the four protons carries the regular tetrahedron into its rotational orbit by a proper rotation
  (the odd relabellings once starred); the 24 rotations close as the cube group O (1 identity, 6 quarter-turns,
  8 third-turns, 9 half-turns), the 12 of the even relabellings as T;
- the spin space is 5 A1 + E + 3 T2 under Td(M) and 5 A + 1E + 2E + 3 T under T; the physical states per species and
  parity (A2 pairs with even parity, A1 with odd); the four isomers of Albert et al.'s Table II, T block, with
  d = 1, 1, 1, 3; the invariant of T x T is unique and of Schmidt rank 3 (9 of the 16 spin states); the monodromy
  groups T/ker Gamma have orders 1, 3, 3, 12; 1E is a homomorphism with 1E(g) = omega;
- the fifteen Cartesian displacements are A1 + E + T1 + 3 T2 and the nine normal to the orbit A1 + E + 2 T2;
- D^J restricted to H for J = 0..6, with the physical states of each J by parity;
- the cell of the identity among the twelve rotations of T, in the bi-invariant metric, is the octahedron
  |x| + |y| + |z| <= 1 in Rodrigues coordinates (3000 random rotations); its vertices are the quarter-turns and the
  face centres of the cube of third-turns at (+-1, +-1, +-1); the face x + y + z = 1 is glued to x + y + z = -1 by
  R -> R R_g^-1 with a third of a twist;
- a third of a turn about the bond to hydrogen 1 returns the tetrahedron with 2, 3, 4 cycled and cycles the axes;
- the packet shares (1 + 8 c3 + 3 c2)/12, (1 - 4 c3 + 3 c2)/12, 3 (1 - c2)/4 against the projector sums.

