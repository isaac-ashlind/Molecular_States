# Final pass progress (temporary; step 9 deletes .final-pass/)

Status: IN PROGRESS
Live check-in trigger: trig_01XfFJ9hVGVvvc8YXyLGJYKf (fires 2026-10-05T17:41Z)
Audit (step 4): three background Agent auditors, launched 16:11Z, reports to .final-pass/audit/{plates,manuscript,compute}.md. If a report is missing after an interruption, re-run that auditor (prompt: read-only, the five headings in CALIBRATION section 8 for its area).

Follow section 11 of CALIBRATION.md at the start of every turn. Tick a step only in the commit that completes it.

- [x] 0. Start: follow the interruption protocol, arm the check-in, set Status to IN PROGRESS
- [x] 1. Read every process doc, README and build.py; fold the necessary rules into source comments (figurestyle.sty,
         primitives.tex headers, check docstrings); delete HANDOFF.md, docs/author-decisions.md, caption-notes.md,
         verification.md, storyline.md, style.md, docs/history/; rewrite README as a short build note
- [x] 2. Remove reference/
- [x] 3. build.py without scaffold (numbering check folded in, reading the manuscript's aux); delete scaffold/;
         fix every comment that points at deleted files
- [ ] 4. Read-only audit (small workflow, record run id): dead definitions and code, duplicated plate machinery and
         inconsistent styles across plates, table formatting, stale comments, arbitrary constants
- [ ] 5. Shared primitives and plate defragmentation, each plate checked with diffsheet before and after
- [ ] 6. compute/, figures/data/ and checks/ cleanup (no unused emissions, no dead code, accurate docstrings)
- [ ] 7. Manuscript LaTeX: preamble, macros, tables (Section 8, Tables A-D, notation page) formatted perfectly
- [ ] 8. Full python3 build.py: every check passes, no collisions, text inside every plate, manuscript clean
- [ ] 9. Delete .final-pass/, commit "Final pass step 9: done", push, cancel the check-in, short final report

Log (one line per step: step, commit, note):
- 0: started 2026-10-05T16:04Z, backstop deleted, check-in armed
- 1: docs folded and deleted (rules in figurestyle.sty, primitives.tex, check docstrings, camera()/kabsch docstrings), figures/preamble.tex removed (unused), README short. Per-plate disclosures for step 5 headers: git show d40f495:docs/caption-notes.md
- 2: reference/ removed (git history keeps it)
- 3: build.py 5 steps (data, checks, plates+pdf/png+collisions, manuscript+numbering check from the manuscript source and aux); scaffold/, previews, proofs and checks/cvd_preview.py removed; plate PDFs reproducible (no date, no trailer ID)
