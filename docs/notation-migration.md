# Notation migration record (approved instructions, 2026-10-04)

A bounded update of the live build: symbols were inventoried by meaning, then replaced role by role with
anchored edits (every replacement matched its full context and an asserted count; no global substitution).
Prose, claims, section order and scope are unchanged except where a symbol's own definition needed a short
local phrase, listed below.

## Old role -> new role, with locations

| Meaning | Was | Now | Where |
|---|---|---|---|
| configuration matrix, reference, family representative | X, X_0, X_0(q) | bold X, X_0, X_0(a) (`\X`) | manuscript SS1-3, 5, 7, 9-11; figs 1, 2, 3, 5, 6, 11 |
| configuration measure | dX | d bold X | manuscript SS2; fig 2 |
| signed right-action matrix | P_s, P_sigma | bold P (`\Pmat{s}`) | manuscript SS1, 10; figs 1, 3, 5 |
| rotation matrices | R (of r), R_h, R_z, R_y | bold R(r), R_h, R_z, R_y (`\Rmat`) | manuscript SS1, 5, 10, 12; figs 1, 5 |
| chart matrices | A_s(q), M_s(q) | bold A_s(a), M_s(a) | manuscript SS10 |
| spatial action operator | U_s, U_g, U_t, U_b, U_tu | sans-serif U (`\Uop{s}`) | manuscript SS7, 10, 12 |
| isotypic projector | P_Gamma | sans-serif P (`\Pproj{\Gamma}`) | manuscript SS7 |
| chart density | script J, J(q,Q) | j_varphi, j_varphi(a,q) (`\jdens`) | manuscript SS10, 12 |
| shape (reference-family) tuple | q, q' | a, a' | manuscript SS10; fig 10 (the zoom's shape point) |
| normal coordinates | Q, Q_a, Q_+, Q_- | q, q_alpha, q_+, q_- | manuscript SS10, 12; figs 6, 10; make_data kappa_rule string |
| normal index and count | a = 1..d | alpha = 1..n_q | manuscript SS10, 12 (R^d, (2 pi hbar)^{d/2}, p in R^d -> n_q) |
| normal-frame vectors | e_a(q), e_b(q') | bold e_alpha(a), e_beta(a') | manuscript SS10; fig 10 |
| umbrella coordinate | eta, eta_0 | iota, iota_0 | manuscript SS10; figs 6, 10 |
| amino offset in the reference family | a in (eta, +-a, z_A) | b_A in (iota, +-b_A, z_A) | manuscript SS10 footnote |
| compensating angle of tu | alpha | omega_tu | manuscript SS10 |
| generic angle in the Wigner identity | alpha | omega | manuscript SS12 |
| lab rotation (group element) | R, T_R, Rr | r_0, T_{r_0}, r_0 r | manuscript SS10 |
| packet width (density SD per coordinate) | lambda (interim; sigma before) | Delta | manuscript SS9; fig 9; verify_symbolic, verify, summary.json key `Delta_over_d` |
| ambient Hilbert space | H_ambient | H_amb (`\hamb`) | manuscript SS2, 3 |
| parity sectors | H^+, H^-, H^pm | H_+, H_-, H_pm | manuscript SS3 |
| free domain of figure 11 | V | calligraphic V | cover, fig 11 |
| spin action in figure 3 | P_sigma Psi(X) | sigma . Psi(bold X), as the text writes it | fig 3 |

Retained on purpose: sigma, sigma* (permutation, permutation-inversion); t, u, b; E, E*; J, M, K; m, kappa; rho;
w_Gamma; chi_stat, chi_pm; A_1, A_2, E upright; delta, ell, theta (water); the type index alpha of SS1 and SS3
(see flags); the macro names in figures/data/numbers.tex (`\etaZero`) and the pic names (`newman-tau-0-eta-minus`),
which are internal paths; the internal Python variable names tau, eta.

Shared macros (the same block in docs/manuscript-reorganized.tex and figures/shared/figurestyle.sty):
`\X`, `\Pmat{s}`, `\Rmat`, `\Uop{s}`, `\Pproj{\Gamma}`, `\jdens`, `\hamb`. Nothing else was added.

## Local phrases adjusted (shown in the diff)

- SS1: "via 3x3 rotation matrices" now names the matrix, "rotation matrices R(r)", so that `r . X = R(r) X` is defined.
- SS10: "Use r for orientation, a for internal shape, and q for normal coordinates" (the letters only).
- SS10: "A lab rotation r_0 gives T_{r_0}(r,a,q) = (r_0 r, a, q)" (the letter only; see flags).

## Coverage

- Prose containing symbols: SS1-3, 5, 7, 9-12 of docs/manuscript-reorganized.tex (the abstract and SS4, 6, 8 carry no
  migrating symbol).
- Bullets: SS7, 9, 10, 11, 12 bullets edited where they name a migrating symbol.
- Display and inline mathematics: 60 anchored replacements, 81 sites; after the pass a search for the stale forms
  (X unbolded in math, P_, U_, R_z, script J, eta, lambda, Q_pm, R^d, d/2, (r,q,Q), dQ, H^pm, "ambient" subscript)
  finds none.
- Captions: the manuscript's captions are placeholders ("author caption; see docs/caption-notes.md"); the caption
  notes were migrated (figures 1, 3, 6, 9, 10, 11) and carry a head note on the ASCII spelling of bold and sans serif.
- Tables: the manuscript has none; the verification table in docs/verification.md was migrated (Delta).
- Figure text: cover (cal V), figure 1 (bold X, R, P_sigma), 2 (bold X, dX), 3 (bold X, P_sigma; sigma . Psi), 5
  (bold X_0, P_b, R_b, R_y), 6 (iota, q, bold X_0(tau, iota)), 9 (Delta/d twice), 10 (iota_0, a, q, iota, bold
  e_pm(a)), 11 (bold X in both columns, the graph and the heads; [bold X]; cal V). Figures 4, 7, 8, 12 carry no
  migrating symbol (figure 12's bold J is the angular-momentum vector).
- Generated assets: all twelve plates, the cover, the previews, the proof sheet and figures/pdf were rebuilt from a
  clean `python3 build.py`; summary.json carries `Delta_over_d` and the migrated kappa rule; the Gaussian data file
  header (`x wA1 wE wA2`) names no symbol and is unchanged.
- Checks: verify_symbolic.py uses the symbol Delta for the width and prints it; verify.py's quadrature check is
  parameterized by Delta. All checks pass (build/build.log).
- Not touched: reference/ (archival snapshots and PDFs), docs/PLAN-v2.md and docs/calibration-record.md (records),
  studies/ (rendered studies, kept as they were pitched).

## Gaussian convention

The live source already used the density convention (position variance lambda^2 of |g|^2, c = exp(-d^2/8 lambda^2),
verified symbolically with amplitude exp(-|x|^2/4 lambda^2)). Delta replaces lambda with the same formula and the
same numbers: no rescaling, the three rows of figure 9 (Delta/d = 1/6, 1/3, 2/3) and the shares (0.341, 0.550,
0.837 for A_1) are unchanged.

## Flags for the author (not decided by the edit)

1. Type index alpha (SS1, SS3: types alpha with multiplicity n_alpha, I_alpha, sigma_alpha) is retained; the normal
   index alpha lives only in SS10 and the two never meet in one derivation. Move the type index to beta if the
   reuse across sections is unwanted.
2. Lab rotation: the instruction's grammar makes group elements lowercase, so the lab rotation is written r_0 and
   its coordinate action T_{r_0}. Any other free letter works; R is now reserved for matrices.
3. The product r R_z(omega) (orientation r followed by a frame turn) keeps r as the group element and bolds the
   turn as a matrix, the same pattern as r A_s(a) in the chart formulas; SS1 states that r acts through its matrix
   R(r). If the author prefers one font throughout such products, the alternative is R(r) R_z(omega).
4. The free domain V of figure 11 and the cover is a region, so the grammar table makes it calligraphic; the text
   never names it. Easy to revert if the author wants a plain V.
5. The dummy index b in the normal-frame expansion (SS10) collided with the generator b; it is now beta.
6. Stabilizer naming (H as the rotational-orbit stabilizer) and the family-versus-thickened-F distinction are
   untouched, as instructed.

## Mathematical corrections (separate from the mechanical diff)

- Projector coefficient: the live source already reads `P_Gamma = (d_Gamma/|G|) sum_g conj(chi_Gamma(g)) U_g`,
  which equals the instruction's chi_Gamma(g^{-1}) form for unitary representations. No change was made; the
  author may prefer the g^{-1} spelling. For the real characters of G_6 and G_12 the two agree.
- Seam and Wigner sign: `D^J_MK(r R_z(omega)) = e^{-iK omega} D^J_MK(r)` with `(r, tau + 2 pi, q) ~ (r R_z(-2 pi
  rho), tau, q)` gives kappa = m + rho K; this is checked exactly in verify_symbolic.py and unchanged.
- Density: j_varphi(a, q) is written without an r argument, as the live source had it (the Haar factor is
  r-independent); no determinant was substituted.

## Build and test commands

    python3 build.py                      # data, checks, plates, previews, proof sheet, collision report
    bash <scratchpad>/finish.sh           # manuscript PDF (21 pages, no undefined macros, no over/underfull boxes)
