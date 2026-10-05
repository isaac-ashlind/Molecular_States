# Final pass progress (temporary; step 9 deletes .final-pass/)

Status: NOT STARTED
Live check-in trigger: trig_016XrphoRBULEYp82BHpvAiN (backstop armed before the pass, fires 2026-10-05T19:02Z; step 0 replaces it)
Workflow runs: none yet

Follow section 11 of CALIBRATION.md at the start of every turn. Tick a step only in the commit that completes it.

- [ ] 0. Start: follow the interruption protocol, arm the check-in, set Status to IN PROGRESS
- [ ] 1. Read every process doc, README and build.py; fold the necessary rules into source comments (figurestyle.sty,
         primitives.tex headers, check docstrings); delete HANDOFF.md, docs/author-decisions.md, caption-notes.md,
         verification.md, storyline.md, style.md, docs/history/; rewrite README as a short build note
- [ ] 2. Remove reference/
- [ ] 3. build.py without scaffold (numbering check folded in, reading the manuscript's aux); delete scaffold/;
         fix every comment that points at deleted files
- [ ] 4. Read-only audit (small workflow, record run id): dead definitions and code, duplicated plate machinery and
         inconsistent styles across plates, table formatting, stale comments, arbitrary constants
- [ ] 5. Shared primitives and plate defragmentation, each plate checked with diffsheet before and after
- [ ] 6. compute/, figures/data/ and checks/ cleanup (no unused emissions, no dead code, accurate docstrings)
- [ ] 7. Manuscript LaTeX: preamble, macros, tables (Section 8, Tables A-D, notation page) formatted perfectly
- [ ] 8. Full python3 build.py: every check passes, no collisions, text inside every plate, manuscript clean
- [ ] 9. Delete .final-pass/, commit "Final pass step 9: done", push, cancel the check-in, short final report

Log (one line per step: step, commit, note):
