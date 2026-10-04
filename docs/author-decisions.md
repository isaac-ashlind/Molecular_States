# Decisions taken, and decisions left to the author

## Taken during the build (author-approved or author-directed)
- Visual identity locked: two ink weights, solid crescent on spheres, generatrix hatching on curved walls, hatch
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
1. Section order. The approved specification numbers the figures 5 Symmetry groups and versions, 6 Feasible
   configuration space, 7 Symmetry species, 8 Nuclear spin weights, 9 Localized states. The Overleaf export has
   sections 5 Symmetry Groups, 6 Methylamine, 7 Versions and Symmetry Species, 8 Localized States, 9 Feasible
   Configuration Space. File names follow the specification. `docs/integration.patch` places each plate by
   content in the export (figure 5 in Methylamine, figures 7 and 8 in Versions and Symmetry Species, figure 9 in
   Localized States, figure 6 in Feasible Configuration Space); LaTeX then numbers them 5, 6, 7, 8, 9 in that order,
   which differs from the file names for 6 to 9. Either reorder the sections to the specification or accept the
   counter's numbering.
2. Final captions, from docs/caption-notes.md.
3. Whether to keep the schematic K2Rb2 positions or replace them with a computed complex geometry.
4. The illustrative component functions of figure 2 and rho = 0.3 of figure 12, if real values are preferred.
