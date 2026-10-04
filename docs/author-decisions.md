# Decisions taken, and decisions left to the author

## Taken during the build (author-approved or author-directed)
- Visual identity locked: two ink weights, two-tone gouache shadow on spheres (the black crescent was replaced at the author's request), generatrix hatching on curved walls, hatch
  only for cut material, flat gouache fills, one oblique projection, labels on leaders or digits inside hydrogens.
- Rods start on the junction curve where the bond leaves the sphere (author: bonds must not enter the nuclei).
- Versions and transformed molecules drawn solid with digits inside the hydrogen discs (author chose option C).
- Lattices: the whole in grey, the retained part in black, the chain in blue, the same on both sides (figure 5).
- Figure 6: zoom shows the cut face itself, hidden versions as faded dashed dots, Newman row keyed to the paths.
- Figure 10: the band of figure 6 is shown being cut and unrolled; the slice at q is the cut face of figure 6.
- Figure 4: no envelopes; dashed boundaries, thin grouping lines, the 16-subgroup lattice of S.
- Figures 8, 12, front page: chosen from alternative studies (docs/illustration-notes.md).
- Water atlas: delta in units of l with the -l/4 column at -0.25 (author correction).
- DECISION items of PLAN-v2 resolved as: matrix inset kept (1); table replaced by the correlation diagram (8);
  slice kept in figure 10 but as the figure 6 cut face, reconstruction molecule kept (10); formula panel dropped
  (12); guide as one scene without frames.

## Open, for the author
1. Section order. The export's sections 5 to 9 (Symmetry Groups, Methylamine, Versions and Symmetry Species,
   Localized States, Feasible Configuration Space) differ from the figure flow of the approved specification.
   As agreed, the integration reorganizes the sections to the figure flow in four parts: `docs/manuscript-reorganized.tex`
   is the export with `\part` headings, the twelve sections of the specification, and the plates placed; the author's
   text is moved verbatim (the Methylamine text continues the Symmetry Groups and Versions section with no heading of
   its own; Nuclear Spin Weights is a new section holding figure 8 and a comment where the author's text goes);
   `docs/integration.patch` is the diff from the export. It compiles with figures 1 to 12 in section order. The
   storyline across the parts is in docs/storyline.md.
2. Final captions, from docs/caption-notes.md.
3. Whether to keep the schematic K2Rb2 positions or replace them with a computed complex geometry.
4. The illustrative component functions of figure 2 and rho = 0.3 of figure 12, if real values are preferred.
