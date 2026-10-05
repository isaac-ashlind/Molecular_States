"""Rigid methane's orientation ball: the twelve rotations of T in the axis-angle ball, the Voronoi cell of the
identity (an octahedron in Rodrigues coordinates), the face gluing, and projected coordinates for the mockup.
Scratch check; if the plate is adopted this moves into make_data.py / verify_methane.py."""
import itertools, math, sys
import numpy as np
rng = np.random.default_rng(7)

# ---- the reference in the C2 frame: hydrogens on the body diagonals, carbon at the origin
L = 1.087
D = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float)/math.sqrt(3)
X0 = (L*D).T                                      # 3 x 4, columns = H1..H4

def fit(X, Y):
    """proper rotation R minimising |R X - Y| (SVD Procrustes); returns R, residual"""
    U, s, Vt = np.linalg.svd(Y @ X.T)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1,1,d]) @ Vt
    return R, np.abs(R @ X - Y).max()

def axis_angle(R):
    ang = math.acos(max(-1, min(1, (np.trace(R)-1)/2)))
    if ang < 1e-9: return np.zeros(3), 0.0
    if abs(ang-math.pi) < 1e-9:
        # axis from the symmetric part
        M = (R + np.eye(3))/2
        k = np.argmax(np.diag(M)); n = M[:,k]/math.sqrt(M[k,k])
        return n, ang
    n = np.array([R[2,1]-R[1,2], R[0,2]-R[2,0], R[1,0]-R[0,1]])/(2*math.sin(ang))
    return n, ang

def rot(n, ang):
    n = np.asarray(n, float); n = n/np.linalg.norm(n)
    K = np.array([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
    return np.eye(3) + math.sin(ang)*K + (1-math.cos(ang))*K@K

# ---- (a) every relabelling is a rotation of X0 (odd ones with the star); classify
def perm_matrix(p):                      # X P_sigma : column j of X P is column p[j] of X  (nucleus p[j] sits where j was)
    P = np.zeros((4,4))
    for j, pj in enumerate(p): P[pj, j] = 1
    return P
def cycle_type(p):
    seen, ct = set(), []
    for i in range(4):
        if i in seen: continue
        c, j = 0, i
        while j not in seen: seen.add(j); j = p[j]; c += 1
        ct.append(c)
    return tuple(sorted(ct, reverse=True))
def parity(p): return sum(1 for i in range(4) for j in range(i) if p[j] > p[i]) % 2

rows, T, O = [], [], []
for p in itertools.permutations(range(4)):
    Y = X0 @ perm_matrix(p)
    star = parity(p) == 1
    if star: Y = -Y
    R, res = fit(X0, Y)
    assert res < 1e-9, (p, res)
    n, ang = axis_angle(R)
    rows.append((p, star, cycle_type(p), round(math.degrees(ang)), n))
    O.append(R)
    if not star: T.append(R)
assert len(T) == 12 and len(O) == 24
by = {}
for p, star, ct, deg, n in rows: by.setdefault((ct, star, deg), []).append((p, n))
for k in sorted(by): print("class", k, "count", len(by[k]))
# the third-turn by +2pi/3 about the bond to H_k is the 3-cycle fixing k
for p, star, ct, deg, n in rows:
    if ct == (3,1) and deg == 120:
        k = [i for i in range(4) if p[i] == i][0]
        toward = np.dot(n, D[k])
        print(" 3-cycle fixing H%d: axis.D[%d] = %+.3f  (nucleus i sits where j was: %s)" % (k+1, k+1, toward, [q+1 for q in p]))
# which quarter-turns: starred 4-cycles, axes x,y,z
qt = [(p, n) for p, star, ct, deg, n in rows if deg == 90]
print("quarter-turns:", [(tuple(q+1 for q in p), np.round(n,3).tolist()) for p, n in qt])

# ---- (b) Voronoi cell of e among T, in Rodrigues coordinates, is the octahedron |x|+|y|+|z| <= 1
def angle_of(R): return axis_angle(R)[1]
def random_rotation():
    q = rng.normal(size=4); q /= np.linalg.norm(q); w, x, y, z = q
    return np.array([[1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
                     [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
                     [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)]])
bad = 0; N = 20000
for _ in range(N):
    R = random_rotation()
    dists = [angle_of(h.T @ R) for h in T]          # d(R, h) = angle(h^-1 R), bi-invariant
    nearest_e = dists[0] <= min(dists[1:]) + 1e-12  # T[0] is the identity (p = (0,1,2,3))
    n, ang = axis_angle(R); rod = math.tan(ang/2)*n
    inside = np.abs(rod).sum() <= 1
    if abs(np.abs(rod).sum() - 1) < 2e-3: continue   # skip the boundary
    if nearest_e != inside: bad += 1
print("Voronoi cell = octahedron |x|+|y|+|z|<=1 (Rodrigues): mismatches", bad, "of", N)
assert bad == 0

# ---- (c) the face x+y+z=1 is the bisector of e and g = R((1,1,1), 2pi/3); right multiplication by g^-1 glues it
#          to the face x+y+z=-1, sending the vertex R_x(pi/2) to R_z(-pi/2): a third of a twist, not straight across
g = rot([1,1,1], 2*math.pi/3)
Rx = rot([1,0,0], math.pi/2); Rz_m = rot([0,0,1], -math.pi/2)
assert np.allclose(Rx @ g.T, Rz_m), "vertex map"
for _ in range(200):
    a, b = rng.random(2); 
    if a + b > 1: a, b = 1-a, 1-b
    rod = a*np.array([1,0,0]) + b*np.array([0,1,0]) + (1-a-b)*np.array([0,0,1])   # on the face x+y+z=1
    ang = 2*math.atan(np.linalg.norm(rod)); R = rot(rod, ang)
    assert abs(angle_of(R) - angle_of(g.T @ R)) < 1e-9                               # equidistant from e and g
    n2, ang2 = axis_angle(R @ g.T); rod2 = math.tan(ang2/2)*n2
    assert abs(rod2.sum() + 1) < 1e-9 and abs(ang2 - ang) < 1e-9                    # lands on x+y+z=-1, same angle
print("face gluing: x+y+z=1 -> x+y+z=-1 by R -> R g^-1, vertex R_x(pi/2) -> R_z(-pi/2)  ok")
# the vertices of the cell are the quarter-turns = starred 4-cycles (checked above: six of them on +-x,+-y,+-z)
# edge midpoint (1/2,1/2,0): angle 2 atan(1/sqrt2)
print("cell edge midpoint angle %.2f deg; vertex (quarter-turn) 90; face centre 60" % math.degrees(2*math.atan(1/math.sqrt(2))))

# ---- (d) projected coordinates for the mockup (orthographic; az, el in degrees)
az, el = math.radians(float(sys.argv[1]) if len(sys.argv)>1 else 18), math.radians(float(sys.argv[2]) if len(sys.argv)>2 else 22)
u = np.array([-math.sin(az), math.cos(az), 0]); v = np.array([-math.sin(el)*math.cos(az), -math.sin(el)*math.sin(az), math.cos(el)])
d = np.array([math.cos(el)*math.cos(az), math.cos(el)*math.sin(az), math.sin(el)])
s = 0.8   # cm per radian
def proj(p): p = np.asarray(p, float); return p@u, p@v, p@d
out = []
out.append("\\def\\ballR{%.4f}\\def\\ballE{%.4f}" % (math.pi*s, math.sin(el)))
hb = [k for k in range(4) if D[k]@d < 0]; hf = [k for k in range(4) if D[k]@d >= 0]
out.append("\\def\\Hback{" + ",".join("%s/%d" % ("abcd"[k], k+1) for k in hb) + "}\\def\\Hfront{" + ",".join("%s/%d" % ("abcd"[k], k+1) for k in hf) + "}")
for k in range(4):
    x, y, z = proj(D[k]*0.66); out.append("\\def\\Hx%s{%.4f}\\def\\Hy%s{%.4f}\\def\\Hd%s{%.4f}" % ("abcd"[k], x, "abcd"[k], y, "abcd"[k], z))
    x, y, z = proj(D[k]*2*math.pi/3*s); out.append("\\def\\Cx%s{%.4f}\\def\\Cy%s{%.4f}\\def\\Cd%s{%.4f}" % ("abcd"[k], x, "abcd"[k], y, "abcd"[k], z))
    x, y, z = proj(-D[k]*2*math.pi/3*s); out.append("\\def\\Dx%s{%.4f}\\def\\Dy%s{%.4f}\\def\\Dd%s{%.4f}" % ("abcd"[k], x, "abcd"[k], y, "abcd"[k], z))
for name, ax in zip("xyz", np.eye(3)):
    for sgn, tag in ((1, "p"), (-1, "m")):
        x, y, z = proj(sgn*ax*math.pi*s); out.append("\\def\\S%s%sx{%.4f}\\def\\S%s%sy{%.4f}\\def\\S%s%sd{%.4f}" % (name, tag, x, name, tag, y, name, tag, z))
        x, y, z = proj(sgn*ax*math.pi/2*s); out.append("\\def\\V%s%sx{%.4f}\\def\\V%s%sy{%.4f}\\def\\V%s%sd{%.4f}" % (name, tag, x, name, tag, y, name, tag, z))
# cell faces: normals (+-1,+-1,+-1); front-facing if n.d > 0; edges visible if on a front face
faces = [np.array(f) for f in itertools.product((1,-1), repeat=3)]
front = [tuple(f) for f in faces if f@d > 0]
def vname(axis, sgn): return "V%s%s" % ("xyz"[axis], "p" if sgn > 0 else "m")
vis, hid = set(), set()
for f in faces:
    verts = [(i, int(f[i])) for i in range(3)]
    for a, b in itertools.combinations(verts, 2):
        e = tuple(sorted([vname(*a), vname(*b)]))
        (vis if tuple(f) in front else hid).add(e)
hid -= vis
out.append("\\def\\cellVisible{" + ",".join("%s/%s" % e for e in sorted(vis)) + "}")
out.append("\\def\\cellHidden{" + ",".join("%s/%s" % e for e in sorted(hid)) + "}")
out.append("\\def\\cellFront{" + ",".join("%s/%s/%s" % tuple(vname(i, f[i]) for i in range(3)) for f in front) + "}")
open(__import__("os").path.join(__import__("os").path.dirname(__file__), "ball-coords.tex"), "w").write("%% emitted by ball_check.py (scratch); az %g, el %g, 0.8 cm per radian\n" % (math.degrees(az), math.degrees(el)) + "\n".join(out) + "\n")
print("front faces", front)
print("wrote ball-coords.tex")
