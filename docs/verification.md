# Verification

Every number and structural claim drawn in the plates is computed in `compute/` and checked three ways.

| Layer | File | What it establishes |
|---|---|---|
| exact and floating arithmetic | `checks/verify.py` | group orders and cosets; [H,B] has 10 subgroups and 17 covers; [H,S] has 36 subgroups and 73 covers; relations btb = t^-1, bt = t^2 b; spin weights and C[G6/H] multiplicities; P_{st} = P_t P_s; the E-pair matrices; torsion species; Gaussian shares and limits; centring to 1e-12; b = R_y(pi) to 1e-15; family actions to 1e-15; min |gX0 - X0| = 2.30 A; emitted data files agree with recomputation; KRb groups |
| symbolic | `checks/verify_symbolic.py` (SymPy) | the four family actions hold identically in (tau, eta); X m = 0 identically; c = exp(-d^2/8 sigma^2) and the shares exactly; the E matrices exactly; the seam relation kappa = m + rho K exactly |
| independent algebra | `checks/verify.g` (GAP) | orders 2, 6, 12, 24, 240; cosets 3, 6, 20, 120; [H,B] 10/17 with orders [2,4,4,4,6,8,12,12,12,24]; three order-12 extensions, any two generate B; [H,S] 36/73; class sizes; spin weights 12/4/8 (even) and 4/12/8 (odd); C[G6/H] = A1 + E; KRb: |S| = 8, 16 subgroups, channel orders 4 |
| text | `checks/collisions.py` | overlapping word boxes per plate (pdftotext -bbox); remaining reports are glyph pairs of one symbol (the ≅ sign, roots in fractions) and the eclipsed digits of figure 6 |

`python3 build.py` runs all of them before compiling; a failing check stops the build. The last full run is in
`build/build.log`. Declared schematic or illustrative: K2Rb2 positions (figure 4), the component functions a, b
(figure 2), rho = 0.3 (figure 12), the tube as a picture of F (figure 6), the sheets (figure 11).
