# Final pass calibration (read fully before touching anything)

Written 2026-10-05 before a context compact, for the author's final consistency, cohesion and defragmentation pass.
This file and `.final-pass/` are TEMPORARY. They are committed only so that an interruption (container reset, context
compact, usage limit, a turn that stops) cannot lose them. The pass's last step deletes `.final-pass/` entirely.
The same text is in the scratchpad at FINAL_PASS_CALIBRATION.md; the repository copy is authoritative.

## 1. The author's brief (their words, condensed)
- "This should be the final consistency, cohesion, and defragmentation pass."
- "Do not stress about layout at this stage, as the prose will change the structure. But formatting must be perfect
  where writing prose will not change it (within figures, tables, etc.)."
- "We must not have any arbitrary choices in the final build."
- "I don't want a handoff doc, just a clean repo with clean latex for me to work from."
- "Everything must be necessary and meaningful. As always, the goal is lean, clean, and canonical."
- "Subagents are allowed as long as they don't lead to fragmentation." Use read-only auditors; apply every change
  myself, centrally, so one hand and one judgment shape the result.
- "Leave nothing but the prose for me." The author writes prose. Everything else must be finished.

## 2. Scope decisions the author made (2026-10-05, answered questions)
1. Process docs: DELETE `HANDOFF.md`, `docs/author-decisions.md`, `docs/caption-notes.md`, `docs/verification.md`,
   `docs/storyline.md`, `docs/style.md`, `docs/history/`. Fold the few durable rules into source comments
   (`figures/shared/figurestyle.sty`, `figures/shared/primitives.tex` headers; what each check verifies into its
   docstring). README becomes a short build note. Read each doc BEFORE deleting it and keep what is necessary.
2. Generated files: KEEP `figures/pdf/`, `figures/png/` (600 dpi, author asked for them), `figures/data/` committed,
   so the LaTeX compiles at once.
3. Tooling: KEEP `compute/` and `checks/`. DROP `scaffold/` (proof sheet goes; the figure-numbering check moves into
   `build.py`, which must not reference scaffold or docs afterwards).
4. `reference/` (8.2 MB, 94 files, the author's materials): REMOVE from the repository (git history keeps it).

## 3. Repository facts
- Repo `/home/user/Molecular_States`, remote `isaac-ashlind/Molecular_States`, work and push on branch `figures-v2`
  only (`git push -u origin figures-v2`). Never open a PR. No model identifiers in commits or files.
- Commit trailer, every commit:
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_01M1roFcbN7NEKLnbHsomMpo`
- State before the pass: HEAD 878a28d, clean tree, `python3 build.py` passes every check (verify.py, antiprism,
  symbolic, spin, methane, GAP), 14 plates with no text collisions, manuscript 25 pages, 0 undefined refs, 0 bad boxes.
- Build: `python3 build.py` (data, checks, plates, previews, numbering check, proof sheet, PDFs + 600 dpi PNGs to
  figures/pdf and figures/png, collisions, manuscript). Single plate: `python3 build.py --only NN`.
- Manuscript only: `bash .final-pass/tools/ms.sh` (pdflatex, bibtex, 2x pdflatex into build/ms, prints pages, bad boxes,
  undefined refs). The manuscript includes plates from `figures/pdf/`, so rebuild plates before judging the PDF.
- Compile one plate to a scratch dir: `cd figures/src && TEXINPUTS=../shared//:../data//: pdflatex
  -interaction=nonstopmode -halt-on-error -output-directory=<dir> figNN-*.tex`.
- Before/after visual check: `python3 .final-pass/tools/diffsheet.py BEFORE_DIR AFTER_DIR OUTPREFIX plate...` (crops every
  changed region at 300 dpi into contact sheets; view them). Text collisions: `python3 checks/collisions.py *.pdf`.
  Text inside the plate: pdftotext -bbox and check x within the plate's bounding box (0..15 cm).

## 4. Pitfalls already hit (do not repeat)
- A `%` comment appended mid-line comments out every node after it on that line (it hid figure 3's chi label).
  Put comments at line ends only. Word counts from collisions.py do not catch a missing label; the diff sheets do.
- Regex over the manuscript: never match from the first `\caption{` to a later label. Anchor on the label and search
  backward (`rindex`) for the nearest `\caption{`. Always assert a single match; check file length; compile.
- Lengthened in-plate lines can run past the plate edge (15 cm). Measure with pdftotext -bbox.
- `blend mode=multiply` for pale hidden points only over white; over a tint it muddies (figure 10's seam point).
- Closing pages: content length moves page breaks. The author said not to stress layout now.

## 5. Writing rules (author's voice and punctuation)
- Voice = the opening pages (Sections 1 and 2): plain declarative sentences, one idea each, "Consider", "Take", "We",
  terms defined at first use, no flourish, simple words, claims no more than shown.
- Prose I write: no colons, semicolons rare, em dashes never. Position words ("at the top", "on the left") instead of
  "Top:" or "Left:". In-plate notes: lowercase, no full stops, no colons.
- Captions were rewritten in this voice (commit 878a28d). Keep them; only fix if a non-prose fact is wrong.
- Headings: a period after the number only ("1. Configuration Space"); unnumbered titles bare. Closing section title
  "Example and Recovery of the Rigid Case" (Title Case like the other headings). Run-in heads without trailing period.

## 6. Notation decisions (all in the text now; keep consistent everywhere)
- Brackets only where they are the standard convention, never arbitrary: [X] the class of X (a point of the
  quotient); [H,B] and [H,S] subgroup-lattice intervals, defined at first use as {G : H <= G <= B}; C[G/H], C[T]
  permutation modules; [M_s(a)]_{beta alpha} matrix entry; real intervals. No invented symbols (L_H^B was withdrawn).
- gamma_s = rearrangement path from X to sX (Section 5's term and letter); in the quotient it is a loop and keeps the
  name gamma_s (no brackets). Loops on figures 11 and 13: gamma_t, gamma_b, gamma_g.
- iota umbrella coordinate (+-iota_0 the versions' values), tau torsion angle, a = (tau, iota) shape, q normal
  coordinates in R^{n_q} (R^9 for methane; never "disc D^9"), D^J Wigner matrix, omega only an angle (roots of unity
  written e^{2 pi i/3}), Rodrigues tan(omega/2) n, eta_0, eta_j, eta_h packets (g is the group element), s, s' in the
  action rule (t is the generator (123)), f_xi, f_zeta the Section 2 components (basis xi, zeta), ell bond scale,
  Delta packet width, Table A classes by relabellings (123), (12)(34), (1234)*, (12)*.
- H_spin includes every nucleus: methylamine (C^2)^{x5} x C^1 x C^3 = 36 A1 + 12 A2 + 24 E, physical weights 36, 12, 24
  (even) and 12, 36, 24 (odd); figure 8 draws the proton factor (32 lines), its caption says nitrogen triples each.
  Methane H_spin = (C^2)^{x4} x C^1.
- Do not re-raise the camera handedness of the methylamine glyphs (the author dismissed that audit; reader-facing logic
  is unaffected). Keep any existing code comment about it accurate and short.

## 7. Plate style rules (from docs/style.md and decisions; fold the necessary ones into source comments)
- Ink #303438 for every line; grey only for greyed-out objects; colour only where it carries meaning (palette in
  figurestyle.sty: red D56860 version points, blue 397BA8 nitrogen, teal 4F9C97 reference family / sphere / rubidium,
  gold E0B84F retained regions / potassium; hues pure).
- Weights: heavy .9 pt, fine .4 pt, thin .3 pt (leaders, axes behind a molecule, closing cube), hatch .22 pt; hatch
  spacing and phase so no line sits on a parallel edge.
- Action line styles at action weight: solid t, dashed u or tu, beaded dots starred. Stealth heads only where a path has
  a direction; a head stops about 3 pt short of the dot it points at, never under it; arrows centred in their gaps.
- Hidden contour = same stroke dashed (`back`); ghost outlines dashed (`ghost`); a hidden point (every edge at it hidden)
  is pale in a ghost outline under the front lines; projections and lifts between levels dotted (`projection`).
  Exception: figure 10's undisplaced ghost is a simple faded solid outline.
- Leaders .3 pt, start on the object (dot edge), stop short of the label, cross nothing where avoidable, no halos.
- Version points are the two-tone \hatom spheres everywhere (radius .08).
- Text under an object centred on it; row and corner titles left-aligned at the margin; labels 9 pt (`lab`), notes 8 pt
  (`sub`); if one label cannot take the plate's size, the plate's labels come down to it (figure 6 is all note size).
- Every plate 15.2 cm wide, included at \linewidth (one print scale).

## 8. Known defragmentation and cleanliness targets (start here, then audit for more)
- Comments pointing to deleted things: compute/make_data.py:183 and checks/verify.py:120 cite docs/verification.md;
  figurestyle.sty:8 cites studies/; fig12 header cites docs/history; scaffold/outline.tex cites caption-notes. After
  deletions, `grep -rn "docs/\|scaffold\|reference/\|HANDOFF\|caption-notes\|author-decisions\|studies/\|history"`
  must return only legitimate hits (docs/manuscript.tex, docs/refs.bib paths in build.py and README).
- build.py: remove scaffold steps (outline numbering, proof sheet, grayscale proof); move the numbering check in
  (labels fig:1..fig:12 must match section numbers; the guide has no numbered label) by reading build/ms aux.
- figures/preamble.tex: keep only if it is needed and meaningful (it was for including plates elsewhere); decide.
- Duplicated plate machinery: sheet parallelograms are defined locally in fig00 (\sheet), fig11 and fig13 (\csheet);
  local style overrides (fig13 thin/teal/tealthin/gold/under, fig05 covlab, fig04 kline/rbline/karrow); several
  ad hoc font sizes (6.5, 7 pt for digits, graph nodes, strip cells, covlab). Unify into primitives where shared, with
  one canonical definition per element type. Every number in a plate either comes from figures/data or is a declared
  layout coordinate; no unexplained magic constants where a shared style should be used.
- compute/make_data.py and figures/data: no unused emissions; every macro in numbers.tex and every pic in
  molecules.tex used by some plate; no dead functions. checks/: no redundant checks, docstrings say what is verified.
- Manuscript preamble: every package and macro used; no dead macros. Tables (Section 8 weights, Tables A-D, notation
  page) formatted perfectly and consistently (rules, alignment, spacing, math style).
- Plate header comments: concise, accurate, no history narration, no references to deleted files.
- Repository: no stray files, .gitignore sensible, README short (what it is, how to build, layout in a few lines).

## 9. Method (each numbered step is a step in `.final-pass/PROGRESS.md`)
1. Read: README, HANDOFF, docs/*.md (all), build.py, figures/shared/*, a sample of plates, compute/, checks/,
   manuscript preamble. Note what must survive in comments.
2. Optionally run read-only auditors (a small workflow, few agents, structured output) for: unused definitions and
   dead code; duplicated plate machinery and inconsistent styles across plates; table formatting; stale comments.
   I apply every change myself.
3. Apply in themed commits: (a) deletions and folding of docs, README; (b) build.py without scaffold; (c) shared
   primitives and plate defragmentation, each plate verified with diffsheet (expect no visual change unless fixing a
   real inconsistency); (d) data/compute/checks cleanup; (e) manuscript LaTeX cleanup and tables.
4. Full `python3 build.py`: all checks pass, no collisions, plates inside their boxes, manuscript compiles with no
   undefined refs or bad boxes. Commit, push. Final message: short, what changed, what was removed, state.

## 10. Autonomy
The author's prompt after the compact is the last input. Ask no questions; make every judgment call yourself under
these rules, and finish the whole pass (deletions, folding, defragmentation, build, commit, push, short final report).
Backup taken before the pass: commit 878a28d (source zip and manuscript PDF delivered to the author); git history also
keeps every deleted file.

## 11. Interruption protocol (follow it at the start of every turn of the pass, including the first)
1. Location and branch. `cd /home/user/Molecular_States`. If the directory is missing, clone the repository there.
   If the checkout is not on `figures-v2`, run `git fetch origin figures-v2 && git checkout figures-v2`. Then
   `git pull --ff-only origin figures-v2`. The pushed branch is the truth.
2. Read `.final-pass/CALIBRATION.md` (this file) and `.final-pass/PROGRESS.md`. If `.final-pass/` does not exist on
   `figures-v2` and the last commit message starts with "Final pass step 9", the pass is DONE: cancel any pending
   check-in triggers listed in that commit message or in `list_triggers` (names starting "Final pass"), reply in one
   line, and stop.
3. `git status`. Uncommitted changes belong to the first unchecked step. Read `git diff`. If they are a coherent part
   of that step, finish the step. If they are broken or unclear, `git stash push -m "final-pass interrupted"` and
   redo the step from its start. Never discard work silently.
4. Safety net. Keep exactly one pending check-in. At the start of the pass and after every completed step, arm a new
   one with send_later (delay about 90 minutes, name "Final pass check-in", message: "Final pass check-in. Follow the
   interruption protocol in .final-pass/CALIBRATION.md on branch figures-v2 and continue the pass from the first
   unchecked step in .final-pass/PROGRESS.md; if the pass is done, cancel this check-in and stop.") and delete the
   previous one with delete_trigger. Record the live trigger id in PROGRESS.md. A check-in that arrives while the
   pass is already done does nothing but confirm and cancel.
5. Workflows. Record every workflow run id and script path in PROGRESS.md as soon as it launches. After an
   interruption, resume it with Workflow({scriptPath, resumeFromRunId}) instead of re-running; read its journal.jsonl
   before assuming results are lost. If the script path is gone (container reset), re-launch the audit.
6. Commit and push at the end of every step, the PROGRESS.md tick in the same commit, message
   "Final pass step N: <what>" plus the trailer. Small steps; never leave a finished step unpushed.
7. Usage limits. If work stops for a usage limit, the pending check-in resumes it after the window resets. Prefer few,
   well-scoped agents; do the editing myself.
