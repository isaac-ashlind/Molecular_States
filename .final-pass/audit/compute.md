# Compute-layer audit (read-only)

Scope: compute/{geometry,groups,make_data}.py, figures/data/*, checks/*. Method: every macro, pic, .dat column and
summary key grepped in figures/src, figures/shared, docs/manuscript.tex and checks (including `\csname` lookups in fig13,
resolved against the emitted lists); AST scan for unused names; a line tracer over make_data + every Python check (lines
never executed); `python3 compute/make_data.py` rewrites figures/data byte-identically; all checks pass (verify.py,
verify_antiprism, verify_symbolic, verify_spin, verify_methane, GAP). compute/ and checks/ line numbers are exact; plate
line numbers are from the working tree at audit time (plates are being edited in step 5).

Two defects come first:
- **geometry.py:126-148 `kabsch_rotation` returns the transpose of the fitted rotation.** Measured: for Y = R0 X it returns
  R0^T to 1e-17 (residual 1.27 for a 0.3 rad turn, 6.97 for a third-turn). It is right only when R0 is symmetric
  (identity, half-turns), which covers both make_data calls. Fix: `R = transpose(U)` (or build H as sum m y x^T).
  Checked with the fix applied: all 5 emitted files stay byte-identical.
- **build.py:67 + verify.g: a failed GAP `Assert` does not stop the build.** A scratch file with `Assert(0, 1 = 2)`
  run as `gap -q -b f.g` (stdin at EOF) prints the error and exits 0. Fix: `gap -q -b --quitonbreak` (measured exit 1).

## 1. Unused emissions

numbers.tex (118 macros; 27 unused, plus unused list fields). molecules.tex: all 48 pics used (fig01:26 grid loop,
fig02:34 samples loop, fig04:20 krb loop, the rest by name), so no item.
- numbers.tex:7-27 — the 21 depth macros `\ball{C,D}{a,b,c,d}d`, `\ballS{xp,xm,yp,ym,zp,zm}d`, `\ballV{...}d`, `\ballFad` are read nowhere (fig13 reads only `x`/`y` via `\csname \p x/y`). Emitted by make_data.py:341-343 (`emit_point`), called at :346-347, :352-354. Fix: emit x and y only.
- numbers.tex:17,19,23 — `\ballSxmx/y`, `\ballSypx/y`, `\ballSzpx/y` (the 3 near-side half-turns) are unread; fig13:52 draws only `\ballSkinBack` = xp,ym,zm. Emitted by make_data.py:352. Fix: emit `ballS` only where `sgn*dot(v,outv) < 0` (the predicate already used at :355).
- numbers.tex:50 `\hasseNodes` — fields 4 (|K/H|) and 5 (name) are bound but never used (primitives.tex:155, fig00:78, fig05:52); labels come from `\hasseLabels`. make_data.py:652-653. Fix: emit `i/x/y`.
- numbers.tex:47 `\krbNodes` — field 5 (name) unused (fig04:35,40 type the labels by hand). make_data.py:629,635. Fix: drop the field.
- numbers.tex:53 `\bigNodes` — field 5 (name) unused (fig05:42,46). make_data.py:752. Fix: drop the field.
- numbers.tex:46 `\chartCells` — field 1 (coset name) unused (fig10:21). make_data.py:523. Fix: drop the field.
- numbers.tex:62 `\torsionSpecies` — fields 2-3 (cos/sin species) unused: fig12:61-70 works the species out again from `m mod 3`. make_data.py:806-807. Fix: have fig12 read `\cs`/`\ss`, which makes the plate data-driven (preferred), or emit only m.
- gaussian-shares.dat column `wA2` (always 0) — fig09:41-42 plot only wA1, wE. Only verify.py:170-172 reads it, to assert that it is 0. make_data.py:775-780. Fix: drop the column and verify's `z`.
- gaussian-shares.dat rows x = 1.01..2.00 (100 of 200 data rows) fall outside fig09's axis (fig09:40 `xmax=1`, `restrict x to domain=0:1`). make_data.py:777. Fix: `range(2, 101)`.
- component-functions.dat — all columns are read (fig02:41-43), so no item.
- summary.json — no plate and not the manuscript read it. Checks read 10 of the 43 keys: spin_weights, local_copies, interval, interval_HS, normal_frame, versions, generic_X_bX_best_rotation_residual, torsion_species (verify.py:156-167), methane_loop, methane_rigid (verify_methane.py:239-243). The other 33 keys are unread. Fix for each line below: delete the key and any computation that exists only to feed it (see section 2).
- make_data.py:293-295 `water_generic`; :304 `water_grid`; :307 `water_samples`.
- make_data.py:312 `methane_X0`; :410-411 `methane_ball`.
- make_data.py:463 `methylamine_frame`; :464 `methylamine_X0`; :465 `methylamine_X0_Xm`; :469 `water_matrices`.
- make_data.py:478 `b_equals_rotation_residual`; :480 `mirror_plane_normal_is_Rb_axis` (the literal `True`, never computed).
- make_data.py:482 `tau_shapes_deg`; :507 `family_actions`, `family_action_residual`; :524 `chart_cells`.
- make_data.py:557 `mla_disp_page_shift`; :567 `generic_X_Xm`; :572 `min_displacement_by_nontrivial_G12_element`.
- make_data.py:576-577 `bt_equals_t2b`, `tb_ne_bt`, `btb_equals_tinv`.
- make_data.py:582-583 `krb_orders`, `krb_any_two_generate`; :597 `krb_core`, `krb_lobes`; :641 `krb_subgroups`, `krb_covers`.
- make_data.py:683 `bond_interval_generators`; :751 `interval_HS_layout`; :765 `coset_action`.
- make_data.py:782 `gaussian_shown`; :797 `sample_components`; :809 `kappa_rule` (a prose string).
- Unread sub-keys of keys that are read: make_data.py:417 `methane_loop.turn_deg`; :528-530 `normal_frame.gram`, `.cfg_res` (neither asserted anywhere); :759 `interval_HS.orders`; :505 the `g` field of `versions` (verify.py:164 ignores it); :808 the species fields of `torsion_species` (verify.py:167 reads only m). Fix: drop them, or assert gram ~ I in verify.py.
- make_data.py:818-820 — prints 11 summary keys to stdout. build.py:55 runs make_data through `run(..., quiet=True)`, which captures and discards the output on success. Fix: print one `wrote` line or nothing.
- Cross-reference (plates): fig01:12, fig03:14, fig06:17, fig07:11, fig08:13, fig09:12 and fig11:16 `\input{numbers.tex}` but use none of its macros. Fix: drop those inputs.

## 2. Dead code

- geometry.py:25 `SPIN` — never read. Fix: use it for the factor of the non-permuted nuclei (see 4), or delete it.
- geometry.py:43-45 `matmul` — never called (verify.py:57-58 `mm` is a copy of it). Fix: delete.
- geometry.py:144-146 — the improper-rotation branch of `kabsch_rotation` is taken in 0 of 73 calls (make_data plus checks). Fix: delete it with the transpose fix.
- geometry.py:221-227 — `phi` is returned in MLA_FRAME and read only by the unread `methylamine_frame`. Fix: keep it local.
- geometry.py:340-342 `draw_order` — never called. Fix: delete.
- groups.py:7 `import math` — unused. Fix: delete.
- groups.py:77-89 `cycle_notation` — never called (its `labels` parameter is unused too). Fix: delete.
- make_data.py:21 `fmt(x, nd=4)` — `nd` is never passed (31 calls), and the `-0.0000` guard assumes nd=4. Fix: drop `nd`.
- make_data.py:28 `label_macros` parameter of `pic_code_rods` — never passed. Fix: drop it.
- make_data.py:63-64 (end-on bond `continue`) and :78-79 (short-arrow `else`) — never executed. Fix: delete both branches.
- make_data.py:105 — the `elements`, `carbon` and `nitrogen` parameters of `pic_code_newman` are never read. Fix: drop them.
- make_data.py:132 (hull early return) and :147, :154-162 (the two-point capsule of `rounded_loop`) — never taken; the hulls have 4 and 3 points. Fix: delete the branch and its docstring sentence.
- make_data.py:180-203, :506-507 `derive_family_actions` — its output goes only to unread keys, and the residual `worst` is never asserted. It repeats verify.py:127-140. Fix: delete it.
- make_data.py:495-497 — `for ch in gname.replace('^2', '2'): pass` is a dead loop. :499-501 rebuilds the word table of :207-208. Fix: have `version_positions` return g (or the labels) and drop :494-501.
- make_data.py:232-242 — `patterns(X)` uses only `len(X)` and is called 4 times for the same 2 patterns; `tangent2` (:241) is never read. Fix: build the two patterns once.
- make_data.py:321 — `CAM_MOL = CAM_BALL` is a second name for one camera. Fix: one name.
- make_data.py:476-478, :480 — `Rb`, `res_b` and the literal `True` exist only for unread keys (a copy of verify.py:122-124). Fix: delete them; keep `bX0` for the pic.
- make_data.py:531-535 — `arrow_gain` and `arrows_for` are never called. Fix: delete; keep `amp`, which :547 uses.
- make_data.py:570-577 — `dmin` and the three relation booleans are unread copies of verify.py:149 and :26-27. Fix: delete.
- make_data.py:581 `p12`, `p34` unused; :582-583 copies verify.py:47-48; :588-597 `krb_name`, `elist`, `krb_core`, `krb_lobes` feed only unread keys. Fix: delete; keep the asserts at :586-587, which back `\krbIdxCore`.
- make_data.py:617 — `order_` is never read. Fix: delete.
- make_data.py:655 — `from itertools import combinations` repeats the import at :8. Fix: use `itertools.combinations`.
- make_data.py:684 — `chainnames` is an identity dict. Fix: test `names[K_] in ('H','G_6','G_{12}','B')`.
- make_data.py:697 — `lev` is never read. :715 (`v = xpos[K]` for a node with no neighbours) never runs. Fix: delete both.
- make_data.py:764-765 — `coset_action_table()` is called only for an unread key, and `cnames` is unused. Fix: delete.
- make_data.py:768-771 `local_copies` — built from a hand-typed `perm_char` and class sizes, and read only by verify.py:157. verify.g:77-80 computes it independently. Fix: delete both, or compute `perm_char` from `Q.cosets(G6, H)`.
- make_data.py:781-782, :797, :809 — `x_show`/`c_show`, `sample_components` and `kappa_rule` exist only for unread keys. Fix: delete.
- verify.py:31 `names`, :37 `names_`, :46 `p12`, `p34`, `Es` — bound and never read. Fix: unpack into `_`.
- verify.py:49 — `Q.close([], 4) | {Q.E(4)}` is redundant because `close` already returns {E}. Fix: `Q.close([], 4)`.
- verify.py:100-103 — tautologies with no code under test. :103 compares literal constants with themselves. Fix: delete.
- verify_antiprism.py:12-13 — `itertools` and `Fraction` are imported and unused. Fix: delete.
- verify_symbolic.py:58 — `name and (...)` is a dead subexpression; :60 `Mtot` is unused. Fix: simplify and delete.
- verify.g:55-56 — `vt` is computed in `nameOf` and never used. Fix: delete.
- verify_spin.py:145 — types 3 and 1, which :136 computed as `sym` and `anti`, and drops the methyl A2 terms. Fix: `A1 = mA1*sym + mA2*anti`, `A2 = mA1*anti + mA2*sym`, `E = mE*(sym+anti)`.

## 3. Naming and docstrings

eta to iota (manuscript :476-477, :730-733, :1112 use iota for the umbrella coordinate; eta_0, eta_j, eta_h are packets there):
- geometry.py:205, 207 — comments "X0(tau, eta)" and "(eta, +-a, z_A)". Fix: iota; the footnote (manuscript:733) writes b_A, not a.
- geometry.py:225, 227 `eta0` (local and dict key); :226-227, 243-244 key `a`. Fix: `iota0`, and `b_A`, since the manuscript's a is the shape (tau, iota).
- geometry.py:231-232, 243-244, 249-253, 275-276, 282-283 — the `eta` parameters, "dX0/deta" and "X0(tau,eta)". Fix: iota.
- make_data.py:17 `ETA0`; :181, 184, 196, 206, 209, 219, 229-231, 238-242, 462, 484-485, 505 (`et`), 508, 560-562. Fix: `IOTA0`, iota.
- verify.py:120, 126, 135 — `G.MLA_FRAME['eta0']` and `eta0`. Its docstring at :10 already says iota. Fix: rename together with geometry.
- verify_symbolic.py:10, 13, 19, 45, 47-50, 58 — the symbol and parameter `eta`, and the printed "X0(tau,eta)". The docstring at :4 already says iota. Fix: iota.
- verify_antiprism.py:1, 10, 47, 51, 54-55 — "(tau, eta) cylinder" and the local `eta`. Fix: iota.
- Cross-reference: fig10:17 plate comment "[-eta0,eta0]".
- Effect, measured on a renamed scratch copy: numbers.tex, molecules.tex and both .dat files stay byte-identical, and no macro name changes. summary.json changes only `methylamine_frame.eta0` to `iota0` (and `a` to `b_A`), keys that nobody reads. The only reader of the Python name is verify.py:120,126. No plate reads it.
- make_data.py:229, 247, 258 — "q=(0,eta0)", "e_a(q)" and "X0(q)" use q for the shape. In the manuscript a is the shape and q the normal coordinates (:665, :693-695). Fix: write a.
- make_data.py:784-786 — "-(b xi + a eta)", "b = -a" and "(xi - eta)" are the old names for f_xi, f_zeta and zeta (manuscript:200-201, :225). Fix: "f_zeta = -f_xi", "(xi - zeta)". The matching names `\sampleA*`/`\sampleB*` (:796, read by fig02:28) and the .dat columns `a`, `b` (:790, read by fig02:41-42) would need the plate updated if renamed.

Inaccurate or stale:
- geometry.py:129-133, :142 — "Jacobi SVD-free", "Newton iteration for the symmetric square root", "exact for small turns", "does not converge for third-turns" and "Y ~ R X needs R = U" are all false. It is Higham's Newton iteration for the polar factor, it converges, and R = U^T. Fix: correct the code and describe it in one sentence.
- verify_methane.py:63-64 — repeats the false "does not converge for third-turns". Fix: "an SVD fit, independent of compute/".
- geometry.py:24 — "Nuclear spin numbers used in the figures": SPIN is unused (see 2 and 4).
- geometry.py:171-172 — "the bisector ... points along -y, so the hydrogens are above": the hydrogens are at +y, and the bisector from the oxygen points along +y. Fix: "+y".
- geometry.py:280 — says it "Returns the unit normal vector", but it returns (normal, mass-orthonormal tangent basis). Fix: the docstring.
- geometry.py:302-305 — "declared schematic in the figure": neither fig04's caption (manuscript:341-344) nor its text says schematic; only the source comment fig04:7 does. It also says "envelopes", but the plate draws `\gline` lines. Fix: "a schematic layout; the plate draws grouping lines".
- groups.py:158-159 — "the number of physical states per spatial level": this holds for the proton factor only, while the manuscript's weights are 3 times larger (nitrogen-14). Fix: say "proton factor", or include the factor (see 4).
- make_data.py:4-6 — "nothing is typed in by hand; checks/verify.py re-derives the same quantities independently". In fact the code types in the Td and T tables, the class data and the isomers (:420-454), `perm_char` (:769), `acts` (:509) and the Gaussian constants (:787). Methane is checked by verify_methane, and several checks reuse the same compute functions (see 5). Fix: rewrite.
- make_data.py:29 — "the v2 molecule primitive" is version history. Fix: name the primitives (`\molatom`, `\molrodj`).
- make_data.py:110 "(declared, not geometry)" and :538-540 "declared": nothing in fig06, fig10 or their captions declares the 20-degree Newman twist or the per-direction arrow scaling. Fix: "a drawing choice, not geometry".
- make_data.py:181-184 — the docstring says "find" and "confirmed", but the function only evaluates hand-typed candidates and asserts nothing. It also cites "Section 10", which states only tu (:731-732); t and b are in Section 12 (:865) and u in neither. This goes with the deletion in section 2.
- make_data.py:216 — "among the 6 candidates": the loop tries 12 (6 tau times 2 iota). Fix: 12.
- make_data.py:298 "figure 1(c)", :481 "figure 6d", :560 "figure 11b", :579 "figure 4b", :584 "figure 4c, nested outlines": no plate and no caption has panel letters, and fig04 draws no nested outlines. Fix: name the plate only and drop "nested outlines".
- make_data.py:486 — "fragment loops ... (figure 4)": `\fragloop*` is used by fig05:24-25. Fix: figure 5.
- make_data.py:763 — "(figures 7, 8)": coset_action, spin_weights and local_copies feed no plate. Fix: delete the heading with the code, or say "for the checks".
- make_data.py:313-318, 355, 383 — "(author)", "the ball view the author chose" and "the old picture" narrate history. Fix: state the camera rule only.
- make_data.py:508 — "under the derived actions": `acts` at :509 is typed in, not derived. Fix: derive the cells from the action formulas, or say "typed".
- make_data.py:532 — "drawn 2.5 x longer", but the gain is 3.0 (and the code is dead). Fix: delete.
- make_data.py:574 — "continuation results" is an opaque label for the unread booleans. Fix: delete with the code.
- verify.py:10, 125 — "(manuscript, Section 10)" for the actions of t, u, b and tu: only tu is in Section 10, t and b are in Section 12, and u is stated nowhere. Fix: cite correctly, or "the reference family's actions".
- verify.py:11, 173 — "summary.json ... agree with recomputation": the summary values are compared with literals (:156-167); only the .dat file is recomputed. Fix: the wording.
- verify.py:75 — a stream-of-consciousness comment ("... use the unnormalised vectors:"). Fix: one line.
- verify.g:2 — "a failed Assert aborts": it aborts the file but exits 0 (see the top of this report). Fix: comment and build flag.
- verify.g:3 — lists "the class sizes of G6", but :47 only prints them. Fix: `Assert(0, SortedList(List(cls,Size)) = [1,2,3])`, or drop them from the header.
- verify_symbolic.py:14 — "(see docs)" points at deleted docs; :64 "approved notation" narrates a decision. Fix: delete both phrases.

## 4. Redundant or missing checks

Redundant (no independent value):
- verify.py:156 checks the same `Q.spin_weights()` output against the same literal as :40, only after a JSON round trip. Fix: keep one.
- verify.py:158 repeats :32-33 (the same `Q.interval_data`, through the JSON). Fix: keep one.
- verify.py:64-84 (E-pair matrices in rationals) is a weaker copy of verify_symbolic.py:81-90. Both type T and Bm by hand rather than building them from `Q.coset_action_table`, which verify.py:38 asserts on its own. Fix: build T and Bm from act_t and act_b in one place, and drop the other copy.
- verify.py:104-116 (quadrature for c) is covered exactly by verify_symbolic.py:65-72. It has value only if SymPy is missing, and SymPy is installed here. Fix: keep one.
- verify.py:86-97 (torsion species) uses the same formula and the same `G6_CHARS` as make_data.py:263-281, and :167 never compares the emitted species. Fix: compare `S['torsion_species']` (or `\torsionSpecies`) with this derivation, or delete.
- verify_symbolic.py:59-62 checks that its own `centre()` centres, which holds by construction. Fix: delete.
- verify_symbolic.py:77-78 asserts (3+6c)/9 = (1+2c)/3, a typed arithmetic identity rather than a derivation. Fix: derive w_A1 from the Gram matrix [[1,c,c],[c,1,c],[c,c,1]] with P_A1 = (1/3) sum_j t^j, and w_A2 = 0 the same way.
- verify_spin.py — 12 A1 + 4 A2 + 8 E is its third derivation (after verify.py:40 and verify.g:73). Its unique facts (methyl 4 A1 + 2 E, amino 3 + 1, the product) appear in neither the manuscript nor any plate text, only in the fig08:8 source comment. Fix: extend it to the full H_spin (tensor with I_3 for nitrogen), which closes the gap below, or delete it.
- verify_antiprism.py — its claim (a trigonal antiprism on the cylinder whose symmetry group has order 12 and is isomorphic to G12) appears nowhere: grep for antiprism, D6, stabiliser and cylinder over the manuscript and the plates finds nothing. Its overlap with stated facts (|G12| = 12, |G12/H| = 6) is already asserted by verify.py:25,29 and verify.g:15,23. Fix: delete it (and build.py:57), or state the fact in the text.

Missing:
- **Nitrogen factor 3: checked nowhere.** `groups.spin_weights` covers the protons only (n_protons=5, spin=0.5), and verify.py:40,156, verify.g:73-74 and verify_spin.py:143 all assert 12/4/8 and 4/12/8. Nothing asserts H_spin = 36 A1 + 12 A2 + 24 E (manuscript:573), the weights 36/12/24 and 12/36/24 (:586-596) or "times three" (:568, fig08 caption :605). Fix: in verify.py, multiply by (2*SPIN['C']+1)(2*SPIN['N']+1) and assert {36,12,24} and {12,36,24}; measured, this passes and the dimension is 96 = 2^5 x 3. In verify.g, multiply spinchar by 3. In verify_spin, use kron(I_3).
- GAP failures are invisible, so add `--quitonbreak` at build.py:67 (see the top of this report).
- Closing section, R_{h1h2} = R_{h2} R_{h1} (:928): unchecked; verify_methane.py:80-81 checks closure only, with `any`. Measured: it holds for all 576 pairs and the reversed order fails. Fix: assert it in verify_methane.
- Closing section, |S| = 48 (:911) and |S/H| = 2 (:932): unchecked. Fix: one assert each in verify_methane.
- Closing section, C[T] = A + 1E + 2E + 3T (:946): unchecked. Fix: decompose the regular character with `TT` in verify_methane.
- Closing section, Table B column "on T" (:1031-1039; A1, A2 to A, E to 1E + 2E, T1, T2 to T): unchecked. Fix: restrict the TD characters to the even elements and decompose with `TT`.
- Closing section, Table A axes (:1022, :1024): "half-turn about x, y, z" for (12)(34) and "half-turn about a cube edge" for (12)* are unchecked. verify_methane.py:83 checks only the multiset of angles; the 3-cycles and 4-cycles are pinned at :161-165. Fix: check the axis of each class.
- Closing section, "the ring carries A when 3 | K, 1E or 2E otherwise" (:1010-1011): unchecked. Fix: decompose e^{-2 pi i K/3} over the C3 of g in verify_methane.
- Closing section, "the reading through h -> R_h exchanges 1E with 2E ... same isomers" (:1058-1059): unchecked. Fix: recompute the isomers with conjugated E characters.
- Covered, for completeness: Section 8 12/4/8 (verify.py:40, verify.g, verify_spin); 32 dimensions with E filling 16 (verify.py:42, verify.g:75); chi_stat = 1 on G6 (verify.py:44); the parity pairing (the weights at -1). Closing: X0 in its orbit, O and T, the nine vibrations, 5 A1 + E + 3 T2, Tables B and C (species, weights, isomers, d, m_nuc, T/ker), the shares, the octahedron cell, gluing, the loop, Schmidt rank 3 and 9/16 (verify_methane).

## 5. make_data versus checks (independence)

Computed the same way twice:
- Family actions: the candidate formulas are typed in three times (make_data.py:188-191, verify.py:128-131, verify_symbolic.py:47-50), and a fourth time as `acts` (make_data.py:509), which builds `\chartCells`. No check compares `acts` with the formulas. Fix: delete the make_data copy and derive the cells from one formula table.
- bX0 = R_y(pi) X0: make_data.py:475-477 is identical to verify.py:122-124.
- min over g in G12 of |gX0 - X0|: make_data.py:571 is character-for-character the same expression as verify.py:149.
- btb = t^-1 and bt = t^2 b: make_data.py:575-577 repeats verify.py:26-27.
- KRb "any two generate": make_data.py:583 is identical to verify.py:48.
- Centring: make_data.py:295, 465 repeat verify.py:121, 143.
- make_data.py:764 and verify.py:37 call the same `Q.coset_action_table`. make_data.py:766 and verify.py:39 call the same `Q.spin_weights`; the only independent checks are verify.g and verify_spin. make_data.py:644 and verify.py:31 call the same `Q.interval_data`; the only independent check is verify.g:26-29.
- Torsion species: make_data.py:263-281 and verify.py:88-96 use the same formula and the same character table.
- gaussian-shares.dat: verify.py:171-172 applies the closed form of make_data.py:779-780 again. It catches a stale file, not a wrong formula; verify_symbolic derives the formula.
- Methane loop: make_data.py:413-416 and verify_methane.py:170-174 run the same algorithm on the same `G.rot_axis` and `G.methane_c2`.
- verify.py:166 ("X and bX are not rotation related") relies on make_data.py:565, which uses the transposed kabsch, so the value is only an upper bound. It equals the true minimum here because the optimum is symmetric: 0.7587146406226174 against the SVD value 0.7587146406226176. Fix: the kabsch fix.

Hard-coded on both sides:
- local_copies: `perm_char` (3, 0, 1) is typed at make_data.py:769, and verify.py:157 compares the result with a literal (1, 0, 1).
- The Td character table is typed at make_data.py:421 and again at verify_methane.py:87-91, so a typo shared by both would pass. Fix: assert row orthogonality in verify_methane, or take the table from GAP (`CharacterTable(SymmetricGroup(4))`) in verify.g.
- entangled_fraction: make_data.py:458 types 16, and verify_methane.py:242 compares with the literal [9, 16] instead of the computed `rank` (:230) times `spin_t['T']`. Fix: `[rank*spin_t['T'], int(chi[(0,1,2,3)])]`.

Structural:
- The checks read summary.json (verify.py:154-167, verify_methane.py:238-243), and no plate reads it. numbers.tex and molecules.tex, which the plates do read, are read by no check, so `\chartCells` (from the typed `acts`), `\versionList` and `\torsionSpecies` are never compared with anything. Fix: either have the checks parse the derived lists from numbers.tex, or emit each such macro from the same value that the checks assert.

make_data versus plates (the same constant typed twice):
- make_data.py:69 `rad_a` repeats primitives.tex:56 (`\radH` .20, `\radX` .36, `\radK` .40, `\radRb` .46). Fix: one source.
- make_data.py:787 (104.5, 14.0) is typed again at fig02:25 (the comb bars computed in TeX) and fig02:44. Fix: emit the bar heights.
- fig09:27, 43 type the shares .341, .550 and .837 (w_A1 at Delta/d = 1/6, 1/3, 2/3). make_data emits only x = 1/3, to an unread key (:782), and no check covers these. Fix: emit the three shares to numbers.tex.
- fig04:35 overrides the x of `\krbNodes` 10 and 12 (it swaps G_K and G_in), and fig05:42 overrides `\bigNodes` 6, 8, 15 and 26. What is drawn is therefore not make_data's layout (:613-622, :742-750), and the indices are fragile. Fix: produce that order and those nudges in make_data.
