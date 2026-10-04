# Molecular figure handoff — 2026-10-03
Use with author's separately supplied complete Overleaf export. This package supplies approved planning and implementation seed, not a finished figure suite. Original prose is untouched.

Read FABLE_PROMPT.md, docs/AUTHORITY.md, PLAN.md, STYLE.md, ISSUES.md. References are handoff inputs, not final repo contents. Current PDF is the current formulation; older paper/slides are contextual.

Seed includes a shared style, separate prose-free outline (comments for figure insertion), build helper, algebra/character/Gaussian checks. No dummy scientific artwork. Fable creates actual standalone sources and renders. No LaTeX compiler available during preparation; seed TeX uncompiled.

Dependencies: Python3; TeX with pdflatex, standalone, TikZ/PGF, PGFPlots, AMS packages, graphicx, geometry, caption and fontenc. Fable tests and reports actual versions. Build seed from its directory: python3 checks/verify.py; python3 build.py. Build helper will refuse missing figures. It never fabricates success.

Final production tree must be lean; remove handoff reference archives and superseded scaffolding. Do not change author prose or supply ghostwritten captions.
