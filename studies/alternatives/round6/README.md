# Round 6: the orientation ball (study D)

One picture for the rigid recovery, drawn after the thread was settled: the orientation
space SO(3) of rigid methane as the axis-angle ball, the twelve rotations of T as red dots
(centre; two tetrahedra of third-turns, through the bonds and through the faces; three
half-turns on the skin, each seen twice), and one cell of the twelve as a wire octahedron
(Albert et al.'s octahedral space; vertices the quarter-turns, faces glued in opposite pairs
by the third-turns). `ball_check.py` verifies the geometry (every relabelling a rotation,
the Voronoi cell an octahedron in Rodrigues coordinates, the face gluing with a third of a
twist) and emits `ball-coords.tex`. Build: `study.sh studies/alternatives/round6/ball.tex`.
If adopted, the computation moves into `compute/make_data.py` and `checks/verify_methane.py`.
