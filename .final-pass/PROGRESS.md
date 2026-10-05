# Final pass progress (temporary; step 9 deletes .final-pass/)

Status: IN PROGRESS
Live check-in trigger: trig_01AK23NB5x2yh4TNrU5rXDTa (fires 2026-10-05T18:57Z)
Audit (step 4): three background Agent auditors, launched 16:11Z, reports to .final-pass/audit/{plates,manuscript,compute}.md. If a report is missing after an interruption, re-run that auditor (prompt: read-only, the five headings in CALIBRATION section 8 for its area).

Follow section 11 of CALIBRATION.md at the start of every turn. Tick a step only in the commit that completes it.

- [x] 0. Start: follow the interruption protocol, arm the check-in, set Status to IN PROGRESS
- [x] 1. Read every process doc, README and build.py; fold the necessary rules into source comments (figurestyle.sty,
         primitives.tex headers, check docstrings); delete HANDOFF.md, docs/author-decisions.md, caption-notes.md,
         verification.md, storyline.md, style.md, docs/history/; rewrite README as a short build note
- [x] 2. Remove reference/
- [x] 3. build.py without scaffold (numbering check folded in, reading the manuscript's aux); delete scaffold/;
         fix every comment that points at deleted files
- [x] 4. Read-only audit (small workflow, record run id): dead definitions and code, duplicated plate machinery and
         inconsistent styles across plates, table formatting, stale comments, arbitrary constants
- [x] 5. Shared primitives and plate defragmentation, each plate checked with diffsheet before and after
- [x] 6. compute/, figures/data/ and checks/ cleanup (no unused emissions, no dead code, accurate docstrings)
- [x] 7. Manuscript LaTeX: preamble, macros, tables (Section 8, Tables A-D, notation page) formatted perfectly
- [x] 8. Full python3 build.py: every check passes, no collisions, text inside every plate, manuscript clean
- [ ] 9. Delete .final-pass/, commit "Final pass step 9: done", push, cancel the check-in, short final report

Log (one line per step: step, commit, note):
- 0: started 2026-10-05T16:04Z, backstop deleted, check-in armed
- 1: docs folded and deleted (rules in figurestyle.sty, primitives.tex, check docstrings, camera()/kabsch docstrings), figures/preamble.tex removed (unused), README short. Per-plate disclosures for step 5 headers: git show d40f495:docs/caption-notes.md
- 2: reference/ removed (git history keeps it)
- 3: build.py 5 steps (data, checks, plates+pdf/png+collisions, manuscript+numbering check from the manuscript source and aux); scaffold/, previews, proofs and checks/cvd_preview.py removed; plate PDFs reproducible (no date, no trailer ID)
- 4: plates.md (23+27+21+44+12 items) and manuscript.md (69 items) in .final-pass/audit; compute.md still running at tick time, needed only for step 6 (re-run the compute auditor if the file is missing then)
- 5a (partial): shared refactor, all 14 plates pixel-identical: dead sty/prim machinery removed; solid spheres and no labels the molecule defaults; thin, hatch, ring, dense, digit, under styles; \vpoint, \vpointhidden, \sheet, \hassecovers; ColDotHidden, ColPatch, ColSphere roles; fig12 species from data. Next 5b: visual unifications (weights, greys, fonts, arrows), then headers. compute.md audit is in.
- 5b (partial): one head for heavy arrows (harrow = tpath head), sarrow/farrow; off-scale weights snapped (.7 covers and mode curves heavy, .35 thin, .45/.65 fig13 fine/heavy, joins fine); greys to ColGreyed/ColFaint; gold tints to ColPatch; \inkdot 1.5pt; digits 6.5pt in every disc; Hasse names dense; \leaderlabel from named coordinates at the dot edge; fig13 head 3pt short; loop and note centring. Checked with diff sheets. Next 5c: tube geometry macros, data-driven numbers (shares, counts, Gaussian, fig05 turn), then 5d headers.
- 5c (partial): \tubegeom and cover view (tube numbers read from the tube; the cover solid takes figure 6 proportions); image style in fig06; make_data emits mlaTurn, spin counts, packet rows and marks, Gaussian centre/width; fig07 amplitudes exact. Next: American spelling pass (author request), then 5d headers.
- 5 spelling: American spelling across manuscript, plates, comments and identifiers (ColGrayed, \nmcenter, \compCenter, center()); refs.bib titles untouched; figure 2 labels f_xi, f_zeta in ink (author)
- 5d: all 14 plate headers rewritten (what the plate shows, what is computed, what is declared, the colors; no colons, no history); stale body comments fixed (V, eta0, author/date narration, drafts); shared headers without colons; fig05 axis thin, fig07 plane arc 14 degrees clear like the key. Step 5 done.
- 6a (partial): kabsch transpose fixed (data identical); GAP --quitonbreak; verify_spin checks the full spin space 36/12/24 and the weights by parity; verify_methane checks Table A rows, R_{h1h2} = R_{h2}R_{h1}, C[T], Table B on-T column, the rings, the 1E/2E exchange. Next 6b: exact rewrite (author: cyclotomic, no numerics) of verify_methane, verify_spin, make_data methane, GAP methane block; then dead code, unused emissions, eta->iota.
- 6b (partial): exact rewrite: verify_methane stdlib exact (integer R_h, Q(omega) class, exact Rodrigues cell, exact vibrations and J ladder), verify_spin exact counts, make_data methane exact (QOmega in groups.py, Chebyshev chi_J), GAP methane block (E(3)); verify.py checks numbers.tex (spin counts, kappa, species, version labels, chart cells, shares); verify_symbolic: iota, shares from projectors, seam for all m; verify_antiprism removed; dead code and unused emissions removed (summary.json 9 keys, no depth macros, trimmed lists, chart cells from the versions, lattice nudges in make_data); analytic tangents; build: SymPy required, no numpy.
- 6: done (last comment fixes; typed tables marked)
- 7: manuscript: caption labels with a period; one table style (plaintable: lettered captions above, small, rules, numbers in math, L columns), Section 8 table is Table A, closing tables B-E with \ref; Notation breakable longtable in the same style with a header row; closing section at body size; run-in heads via \paragraph in Title Case with \ref numbers; section labels on all twelve; ties before \cite; \act, \cong, C_2, 2pi/3, pi/3, gamma_s, u_i (no clash with zeta, e); bookmarks and PDF title; title Their; fig05 bounding box trimmed (float fit); 26 pages, 0 bad boxes, 0 undefined
- 8: full build passes: all checks (verify, spin, methane, symbolic, GAP), 14 plates no collisions and no text outside, manuscript 26 pages, 0 undefined, 0 bad boxes, numbering check counts figure captions only
