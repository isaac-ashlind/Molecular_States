"""Generate the verified figure data in figures/data/ from compute/.

Run from the repository root:  python3 compute/make_data.py
Everything emitted here is recomputed from declared masses, model parameters
and group definitions; nothing is typed in by hand.  checks/verify.py
re-derives the same quantities independently and asserts them.
"""
import itertools, json, math, cmath, os, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(__file__))
import geometry as G
import groups as Q

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'figures', 'data')
os.makedirs(DATA, exist_ok=True)

ETA0 = G.MLA_FRAME['eta0']
CAM_MLA = G.camera(azimuth_deg=50.0, elevation_deg=25.0)
CAM_FLAT = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))   # page = xy-plane

def fmt(x, nd=4):
    s = f'{x:.{nd}f}'
    if s.startswith('-0.0000'): s = '0.0000'
    return s

# --------------------------------------------------------------- pic emit ---

def pic_code_rods(X, elements, bonds, cam, label_macros=None, arrows=None, labels=False, only=None):
    """TikZ pic code for the v2 molecule primitive with correct layering.

    Items (atoms and rods) are painted back to front.  A rod is painted just after its farther atom
    and before its nearer atom.  It starts where the bond cylinder leaves the far sphere (the junction
    circle, projected as a half-ellipse computed in TeX from the atom radius) and vanishes under the
    silhouette of the near sphere.  Rod depth = depth of the far atom + a small epsilon.
    `only`, a set of atom indices, restricts the drawing to those atoms and the rods that touch them.
    """
    right, up, out = cam
    P = G.project(X, cam)
    n = len(X)
    if label_macros is None:
        label_macros = ['\\mol' + 'ABCDEFG'[i] for i in range(n)]
    depth = [P[i][2] for i in range(n)]
    rods = []
    for a, b in bonds:
        i, j = a-1, b-1
        far, near = (i, j) if P[i][2] <= P[j][2] else (j, i)
        rods.append((P[far][2] + 1e-6, 1, i, j))
        # a coplanar pair ties in depth: the near sphere must still be painted after the rod
        depth[near] = max(depth[near], P[far][2] + 2e-6)
    items = [(depth[i], 0, i, None) for i in range(n)] + rods
    items.sort(key=lambda t: (t[0], t[1]))
    lines = []
    for z, kind, i, j in items:
        if only is not None and i not in only and (kind == 0 or j not in only):
            continue
        if kind == 0:
            lines.append(f'\\molatom{{{elements[i]}}}{{{label_macros[i]}}}{{{fmt(P[i][0])}}}{{{fmt(P[i][1])}}}')
        else:
            far, near = (i, j) if P[i][2] <= P[j][2] else (j, i)
            u = G.sub(X[near], X[far]); L = math.sqrt(G.dot(u, u)); u = G.scale(u, 1/L)
            ux, uy, uz = G.dot(u, right), G.dot(u, up), G.dot(u, out)
            dn = math.hypot(ux, uy)
            if dn < 1e-6:
                continue            # end-on bond: hidden under the near sphere
            rad = RAD_MACRO.get(elements[far], '\\radX')
            lines.append(f'\\molrodj{{{fmt(P[far][0])}}}{{{fmt(P[far][1])}}}{{{fmt(P[near][0])}}}{{{fmt(P[near][1])}}}'
                         f'{{{fmt(ux/dn)}}}{{{fmt(uy/dn)}}}{{{fmt(uz)}}}{{{rad}}}')
    if arrows is not None:
        rad_a = {'H': 0.20, 'C': 0.36, 'N': 0.36, 'O': 0.36, 'K': 0.40, 'Rb': 0.46}
        for k in range(n):
            if arrows[k] is not None:
                dx, dy = arrows[k]
                # the shaft starts on the atom's silhouette, not at its centre (visible on a light disc)
                L = math.hypot(dx, dy); r = rad_a.get(elements[k], 0.36)
                if L > r + 0.15:
                    ux, uy = dx/L, dy/L
                    lines.append(f'\\molarrow{{{fmt(P[k][0] + r*ux)}}}{{{fmt(P[k][1] + r*uy)}}}{{{fmt(dx - r*ux)}}}{{{fmt(dy - r*uy)}}}')
                else:
                    lines.append(f'\\molarrow{{{fmt(P[k][0])}}}{{{fmt(P[k][1])}}}{{{fmt(dx)}}}{{{fmt(dy)}}}')
    if labels:
        # label direction: away from the mean page position of the bonded neighbours (or straight up)
        for k in range(n):
            nb = [b-1 for a, b in bonds if a-1 == k] + [a-1 for a, b in bonds if b-1 == k]
            if nb:
                mx = sum(P[m][0] for m in nb)/len(nb); my = sum(P[m][1] for m in nb)/len(nb)
                ax, ay = P[k][0]-mx, P[k][1]-my
                nrm = math.hypot(ax, ay)
                ax, ay = (ax/nrm, ay/nrm) if nrm > 1e-6 else (0.0, 1.0)
            else:
                ax, ay = 0.0, 1.0
            lines.append(f'\\mollabel{{{elements[k]}}}{{{label_macros[k]}}}{{{fmt(P[k][0])}}}{{{fmt(P[k][1])}}}{{{fmt(ax)}}}{{{fmt(ay)}}}')
    return lines

RAD_MACRO = {'H': '\\radH', 'C': '\\radX', 'N': '\\radX', 'O': '\\radX', 'K': '\\radK', 'Rb': '\\radRb'}

def pic_code_grouped(X, elements, lines_, cam):
    """Spheres with thin grouping lines between the grouped nuclei (lines first, then atoms, then labels)."""
    P = G.project(X, cam)
    out = []
    for a, b in lines_:
        out.append(f'\\gline{{{fmt(P[a-1][0])}}}{{{fmt(P[a-1][1])}}}{{{fmt(P[b-1][0])}}}{{{fmt(P[b-1][1])}}}')
    out += pic_code_rods(X, elements, [], cam, labels=True)
    return out

def pic_code_newman(X, elements, methyl=(1, 2, 3), amino=(4, 5), carbon=6, nitrogen=7, twist_deg=0.0):
    """Newman projection down the C-N axis (viewer on the N side): back carbon as a circle with its
    hydrogens on fine bonds from the rim; front nitrogen as a disc at the centre with its hydrogens on
    rod bonds from the centre.  Page coordinates are the molecular (x, y) in angstrom; the digit is the
    column label (version labels are supplied by the \\molA.. macros).  twist_deg turns the back set on
    the page by a small angle, the drawing convention for an eclipsed projection (declared, not geometry)."""
    lines = []
    c, s_ = math.cos(math.radians(twist_deg)), math.sin(math.radians(twist_deg))
    for k in methyl:
        x, y, _ = X[k-1]
        x, y = c*x - s_*y, s_*x + c*y
        lines.append(f'\\nmback{{{fmt(x)}}}{{{fmt(y)}}}{{\\mol{"ABCDEFG"[k-1]}}}')
    lines.append('\\nmcircle')
    for k in amino:
        x, y, _ = X[k-1]
        lines.append(f'\\nmfront{{{fmt(x)}}}{{{fmt(y)}}}{{\\mol{"ABCDEFG"[k-1]}}}')
    lines.append('\\nmcentre')
    return lines

def emit_pic(out, name, lines):
    out.append(f'\\tikzset{{pics/{name}/.style={{code={{%')
    for l in lines:
        out.append('  ' + l + '%')
    out.append('}}}')

def hull(points):
    pts = sorted(set(points))
    if len(pts) <= 2: return pts
    def crossp(o, a, b): return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lower = []
    for p in pts:
        while len(lower) >= 2 and crossp(lower[-2], lower[-1], p) <= 0: lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and crossp(upper[-2], upper[-1], p) <= 0: upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def rounded_loop(points, margin):
    """Points of a closed curve around the convex hull, offset outward by margin.

    Two points give a capsule (stadium) sampled at 10 points.

    Each hull vertex is replaced by three points on the circle of radius
    margin about it (along the two edge normals and the bisector), so that
    plot[smooth cycle] through them hugs the hull with rounded corners.
    """
    h = hull(points)
    if len(h) == 2:
        (ax, ay), (bx, by) = h
        dx, dy = bx-ax, by-ay; L = math.hypot(dx, dy); ux, uy = dx/L, dy/L
        out = []
        for ang in (90, 135, 180, 225, 270):
            th = math.radians(ang); out.append((ax + margin*(ux*math.cos(th) - uy*math.sin(th)), ay + margin*(uy*math.cos(th) + ux*math.sin(th))))
        for ang in (270, 315, 0, 45, 90):
            th = math.radians(ang); out.append((bx + margin*(ux*math.cos(th) - uy*math.sin(th)), by + margin*(uy*math.cos(th) + ux*math.sin(th))))
        return out
    m = len(h)
    out = []
    for k in range(m):
        p = h[k]; prev = h[k-1]; nxt = h[(k+1) % m]
        def outward_normal(a, b):
            dx, dy = b[0]-a[0], b[1]-a[1]
            L = math.hypot(dx, dy)
            return (dy/L, -dx/L)   # right-hand normal for a counter-clockwise hull
        n1 = outward_normal(prev, p); n2 = outward_normal(p, nxt)
        bis = (n1[0]+n2[0], n1[1]+n2[1]); bl = math.hypot(*bis)
        bis = (bis[0]/bl, bis[1]/bl) if bl > 1e-9 else n1
        for nn in (n1, bis, n2):
            out.append((p[0]+margin*nn[0], p[1]+margin*nn[1]))
    return out

# -------------------------------------------------------------- derivations --

def derive_family_actions():
    """For g in {t, u, b, tu}: find (tau', eta') and A_g with g.X0(tau,eta) = A_g X0(tau',eta').

    The candidate formulas are those of the manuscript's Section 10; here they
    are confirmed numerically at several (tau, eta) by a residual check.
    """
    results = {}
    cands = {
        't':  (lambda ta, et: (ta + 2*math.pi/3, et), G.identity(3), 'I'),
        'u':  (lambda ta, et: (ta + math.pi, -et), G.rot_z(math.pi), 'R_z(\\pi)'),
        'b':  (lambda ta, et: (-ta, et), G.rot_y(math.pi), 'R_y(\\pi)'),
        'tu': (lambda ta, et: (ta - math.pi/3, -et), G.rot_z(math.pi), 'R_z(\\pi)'),
    }
    els = {'t': Q.t, 'u': Q.u, 'b': Q.b, 'tu': Q.tu}
    worst = 0.0
    for name, (f, A, Aname) in cands.items():
        for ta, et in [(0.0, ETA0), (0.37, 0.6*ETA0), (-1.1, -0.2), (2.0, 0.0)]:
            X = G.methylamine(ta, et)
            gX = G.apply_perm_inversion(X, els[name][0], els[name][1])
            tp, ep = f(ta, et)
            Y = G.rotate(A, G.methylamine(tp, ep))
            worst = max(worst, G.config_distance(gX, Y))
        results[name] = Aname
    return results, worst

def version_positions():
    """(tau_g, eta_g) for the six coset representatives of G12/H."""
    reps = [('H', [], 'E'), ('tH', ['t'], 't'), ('t^2H', ['t', 't'], 't^2'),
            ('uH', ['u'], 'u'), ('tuH', ['t', 'u'], 'tu'), ('t^2uH', ['t', 't', 'u'], 't^2u')]
    X0 = G.methylamine(0.0, ETA0)
    out = []
    for cname, word, gname in reps:
        g = Q.E(Q.N)
        for w in word:
            g = Q.mul(g, {'t': Q.t, 'u': Q.u}[w])
        gX = G.apply_perm_inversion(X0, g[0], g[1])
        # search the derived position among the 6 candidates
        best = None
        for ta in [k*math.pi/3 for k in range(6)]:
            for et in (ETA0, -ETA0):
                R, res = G.kabsch_rotation(G.methylamine(ta, et), gX, G.MLA_MASSES)
                if best is None or res < best[0]:
                    best = (res, ta, et)
        res, ta, et = best
        assert res < 1e-9, (cname, res)
        out.append((cname, gname, ta, et))
    return out

def normal_frame_data():
    """Two symmetry-adapted normal vectors at q=(0,eta0) and their images under tu."""
    tau, eta = 0.0, ETA0
    tau2, eta2 = tau - math.pi/3, -eta
    def patterns(X):
        n = len(X)
        zero = (0.0, 0.0, 0.0)
        stretch = [zero]*n; stretch[5] = (0, 0, -1.0); stretch[6] = (0, 0, 1.0)
        twist = [zero]*n; twist[3] = (0, 0, 1.0); twist[4] = (0, 0, -1.0)
        return stretch, twist
    X = G.methylamine(tau, eta); Xp = G.methylamine(tau2, eta2)
    e_plus, tangent = G.normal_vector(patterns(X)[0], tau, eta)
    e_minus, _ = G.normal_vector(patterns(X)[1], tau, eta)
    ep_plus, tangent2 = G.normal_vector(patterns(Xp)[0], tau2, eta2)
    ep_minus, _ = G.normal_vector(patterns(Xp)[1], tau2, eta2)
    m = G.MLA_MASSES
    # orthonormality and mass-orthogonality to tangents
    gram = [[G.mass_inner(a, b, m) for b in (e_plus, e_minus)] for a in (e_plus, e_minus)]
    ortho = max(abs(G.mass_inner(e, tv, m)) for e in (e_plus, e_minus) for tv in tangent)
    # transformation under tu: e_a(q) P_tu = A sum_b e_b(q') M_ba, with A = R_z(pi)
    A = G.rot_z(math.pi); Ainv = G.transpose(A)
    M = [[0.0, 0.0], [0.0, 0.0]]; resid = 0.0
    for a, ea in enumerate((e_plus, e_minus)):
        img = G.apply_perm_inversion(ea, Q.tu[0], Q.tu[1])
        img = G.rotate(Ainv, img)
        recon = [(0.0, 0.0, 0.0)]*len(img)
        for bb, eb in enumerate((ep_plus, ep_minus)):
            M[bb][a] = G.mass_inner(eb, img, m)
            recon = [G.add(r, G.scale(v, M[bb][a])) for r, v in zip(recon, eb)]
        resid = max(resid, G.config_distance(img, recon, m))
    # also confirm the configuration part: X0(q) P_tu = A X0(q')
    cfg_res = G.config_distance(G.apply_perm_inversion(X, Q.tu[0], Q.tu[1]), G.rotate(A, Xp))
    return dict(e_plus=e_plus, e_minus=e_minus, gram=gram, ortho=ortho, M=M,
                resid=resid, cfg_res=cfg_res, tangent_dim=len(tangent))

def torsion_mode_species():
    """Species of cos(m tau), sin(m tau) under G6 with t: tau -> tau + 2pi/3, b: tau -> -tau."""
    out = []
    for m in range(0, 7):
        if m == 0:
            out.append((0, 'A1', None)); continue
        # rep on (cos, sin): U_t = rotation by -2 pi m/3, U_b = diag(1,-1)
        chi_t = 2*math.cos(2*math.pi*m/3); chi_b = 0.0; chi_E = 2
        mult = {}
        for name, chi in Q.G6_CHARS.items():
            mult[name] = (chi['E']*chi_E + 2*chi['t']*chi_t + 3*chi['b']*chi_b)/6
        mult = {k: int(round(v)) for k, v in mult.items()}
        if mult == {'A1': 0, 'A2': 0, 'E': 1}:
            out.append((m, 'E', 'E'))
        elif mult == {'A1': 1, 'A2': 1, 'E': 0}:
            out.append((m, 'A1', 'A2'))     # cos even under b -> A1, sin odd -> A2
        else:
            raise AssertionError((m, mult))
    return out

# ------------------------------------------------------------------ main ---

def main():
    summary = {}
    mol = ['% Generated by compute/make_data.py. Coordinates in angstrom; depth-sorted.',
           '% Each \\molatom carries a label macro so figures can relabel versions.']
    num = ['% Generated by compute/make_data.py. Numeric labels and lists for the figures.']

    # ---- water
    wX = G.water(0.20, 112.0)          # a generic, visibly asymmetric configuration
    summary['water_generic'] = {'delta_ratio': 0.20, 'theta_deg': 112.0,
                                'X': [[round(v, 4) for v in c] for c in wX],
                                'Xm': [round(v, 12) for v in G.mass_moment(wX, G.WATER_MASSES)]}
    bonds_w = [(1, 3), (2, 3)]
    Rw = G.rot_z(math.radians(40.0))
    # shape grid for figure 1(c): delta ratios and angles, drawn mass-centred
    grid_d = [-0.5, -0.25, 0.0, 0.25, 0.5]; grid_t = [60.0, 90.0, 120.0, 150.0, 180.0]
    for i, d in enumerate(grid_d):
        for j, th in enumerate(grid_t):
            Xg = G.water(d, th)
            emit_pic(mol, f'mol3d-water-grid-{i}{j}', pic_code_rods(Xg, G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    summary['water_grid'] = {'delta_ratios': grid_d, 'theta_deg': grid_t}
    # three sample configurations for figure 2
    samples = [(0.0, 90.0), (0.0, 110.0), (0.0, 140.0)]    # on the delta = 0 slice, in order of theta
    summary['water_samples'] = samples

    # ---- methane (the closing pages): the rigid tetrahedron in its C2 frame, the orientation ball SO(3) with the twelve
    # symmetry rotations and one cell, the vibrational species, the J ladder; the numbers against Albert et al.
    CH4 = G.methane_c2()
    summary['methane_X0'] = [[round(v, 4) for v in c] for c in CH4]
    # one camera for both panels, so that bond 1 on the molecule points along the arrow in the ball (author). It is the
    # suite's G.camera(108, 21), the ball view the author chose, with the x and z components of its three vectors
    # exchanged. G.camera builds a left-handed frame (a mirror image); the exchange x <-> z is a mirror symmetry of the
    # cube, the dual octahedron and the bond tetrahedron (it fixes bond 1 and exchanges hydrogens 2 and 4), so the
    # frame becomes right-handed, every ball element lands where the old picture drew one of its kind, the arrow lands
    # on the third-turn about bond 1, and the molecule is a true view with x to the right, y up, z toward the viewer.
    def _swap_xz(v): return (v[2], v[1], v[0])
    CAM_BALL = tuple(_swap_xz(v) for v in G.camera(azimuth_deg=108.0, elevation_deg=21.0))
    CAM_MOL = CAM_BALL
    assert G.dot(G.cross(CAM_BALL[0], CAM_BALL[1]), CAM_BALL[2]) > 0.999          # right-handed
    assert _swap_xz(tuple(CH4[0])) == tuple(CH4[0]) and _swap_xz(tuple(CH4[1])) == tuple(CH4[3])   # x <-> z fixes H1 and exchanges H2 and H4
    emit_pic(mol, 'mol3d-ch4-c2', pic_code_rods(CH4, G.CH4_ELEMENTS, G.CH4_BONDS, CAM_MOL, labels=True))
    front, back = [], []
    for name, v in (('X', (1.0, 0.0, 0.0)), ('Y', (0.0, 1.0, 0.0)), ('Z', (0.0, 0.0, 1.0))):
        px, py, pz = G.project([v], CAM_MOL)[0]
        num.append(f'\\def\\chAxis{name}x{{{px:.4f}}}\\def\\chAxis{name}y{{{py:.4f}}}')
        (front if pz >= 0 else back).append(name + '/1'); (back if pz >= 0 else front).append(name + '/-1')   # each half-axis: toward the viewer over the molecule, away from it behind
    num.append('\\def\\chAxisFront{' + ','.join(front) + '}\\def\\chAxisBack{' + ','.join(back) + '}')
    # the orientation space as the axis-angle ball (direction the axis, distance the angle, radius pi; the skin is the
    # half-turns, each point its own antipode). The twelve rotations of T: the centre, the eight third-turns at 2pi/3
    # along the body diagonals, the three half-turns on the skin along the C2 axes. The cell of the centre is drawn as
    # the schematic exact in Rodrigues coordinates: the octahedron dual to the cube spanned by the eight third-turns,
    # its vertices (the quarter-turns) on the cube's face centres, every edge straight; in the ball the true edges bow
    # outward (the caption says so). The cube itself is a reading aid: no edge between two points is intrinsic.
    BALL_R = 2.0                                                # cm for the angle pi: the ball's radius on the page
    bs = BALL_R / math.pi                                       # cm per radian
    D = [G.unit(CH4[k]) for k in range(4)]                      # bond directions = body diagonals
    outv = CAM_BALL[2]
    def emit_point(name, v, scale=1.0):
        px, py, pz = G.project([G.scale(v, scale)], CAM_BALL)[0]
        num.append(f'\\def\\{name}x{{{px:.4f}}}\\def\\{name}y{{{py:.4f}}}\\def\\{name}d{{{pz:.4f}}}')
    num.append(f'\\def\\ballR{{{BALL_R:.4f}}}')
    for k, tag in enumerate('abcd'):
        emit_point('ballC' + tag, D[k], 2 * math.pi / 3 * bs)   # third-turn about bond k, +2pi/3: a cube vertex
        emit_point('ballD' + tag, D[k], -2 * math.pi / 3 * bs)  # the opposite turn
    axes = {'x': (1.0, 0.0, 0.0), 'y': (0.0, 1.0, 0.0), 'z': (0.0, 0.0, 1.0)}
    face_centre = 2 * math.pi / 3 * bs / math.sqrt(3.0)          # the cube's half-edge: where its face centres sit on the axes
    for a_, v in axes.items():
        for sgn, t in ((1.0, 'p'), (-1.0, 'm')):
            emit_point(f'ballS{a_}{t}', v, sgn * math.pi * bs)         # half-turn, on the skin
            emit_point(f'ballV{a_}{t}', v, sgn * face_centre)           # quarter-turn, drawn on the cube's face centre: a vertex of the dual octahedron
    emit_point('ballFa', D[0], 2 * math.pi / 9 * bs)             # where the lift leaves the octahedron: the centroid of its face toward bond 1 (the dual octahedron's face plane x+y+z = c lies at c/sqrt3 on the diagonal)
    num.append('\\def\\ballSkinBack{' + ','.join(f'{a_}{t}' for a_, v in axes.items() for sgn, t in ((1.0, 'p'), (-1.0, 'm')) if sgn * G.dot(v, outv) < 0) + '}')   # the far representatives, drawn pale (author: so they look far away)
    # the octahedron: faces with normals (+-1,+-1,+-1); an edge is visible if it lies on a face turned to the viewer
    faces = [f for f in itertools.product((1, -1), repeat=3)]
    front = [f for f in faces if G.dot(f, outv) > 0]
    def vname(i, sg): return 'ballV' + 'xyz'[i] + ('p' if sg > 0 else 'm')
    vis, hid = set(), set()
    for f in faces:
        vs = [vname(i, f[i]) for i in range(3)]
        for e in itertools.combinations(vs, 2):
            (vis if f in front else hid).add(tuple(sorted(e)))
    hid -= vis
    num.append('\\def\\ballCellVisible{' + ','.join(f'{a1}/{b1}' for a1, b1 in sorted(vis)) + '}')
    num.append('\\def\\ballCellHidden{' + ','.join(f'{a1}/{b1}' for a1, b1 in sorted(hid)) + '}')
    # the cube: edges join third-turns differing in one sign; visible if on a cube face turned to the viewer
    verts = {('ballC' + 'abcd'[k]): tuple(D[k]) for k in range(4)}
    verts.update({('ballD' + 'abcd'[k]): tuple(-c for c in D[k]) for k in range(4)})
    sgn = {n_: tuple(1 if c > 0 else -1 for c in v) for n_, v in verts.items()}
    cvis, chid = set(), set()
    for a1, b1 in itertools.combinations(sorted(verts), 2):
        diff = [i for i in range(3) if sgn[a1][i] != sgn[b1][i]]
        if len(diff) != 1: continue
        shared = [i for i in range(3) if i not in diff]
        visible = any(sgn[a1][i] * outv[i] > 0 for i in shared)
        (cvis if visible else chid).add((a1, b1))
    num.append('\\def\\ballCubeVisible{' + ','.join(f'{a1}/{b1}' for a1, b1 in sorted(cvis)) + '}')
    num.append('\\def\\ballCubeHidden{' + ','.join(f'{a1}/{b1}' for a1, b1 in sorted(chid)) + '}')
    assert len(cvis) + len(chid) == 12
    # a point is hidden when every edge at it is hidden: one corner of the cube, one vertex of the cell; drawn pale,
    # under the front lines, like the far half-turns (author)
    def split_points(names, vis_e, hid_e):
        back = [n_ for n_ in sorted(names) if all(e in hid_e for e in vis_e | hid_e if n_ in e)]
        return [n_ for n_ in sorted(names) if n_ not in back], back
    cube_front, cube_back = split_points(verts, cvis, chid)
    cell_front, cell_back = split_points({vname(i, sg) for i in range(3) for sg in (1, -1)}, vis, hid)
    assert len(cube_back) == 1 and len(cell_back) == 1, (cube_back, cell_back)
    assert all(G.dot(verts[n_], outv) < 0 for n_ in cube_back)
    for tag, lst in (('CubeFront', cube_front), ('CubeBack', cube_back), ('CellFront', cell_front), ('CellBack', cell_back)):
        num.append(f'\\def\\ball{tag}Pts{{' + ','.join(lst) + '}')
    # the three great circles of the coordinate planes as the depth cue: front and back arcs
    def circle_pts(plane):
        out = []
        for t in range(0, 720, 4):
            cth, sth = math.cos(math.radians(t)), math.sin(math.radians(t))
            v = {'XY': (cth, sth, 0.0), 'YZ': (0.0, cth, sth), 'ZX': (sth, 0.0, cth)}[plane]
            out.append(G.project([tuple(BALL_R * x for x in v)], CAM_BALL)[0])
        return out
    for plane in ('XY', 'YZ', 'ZX'):
        pts = circle_pts(plane)
        for tag, sel in (('Front', lambda z: z >= 0), ('Back', lambda z: z < 0)):
            start = next(i for i in range(len(pts)) if sel(pts[i][2]) and not sel(pts[i - 1][2]))
            run = []
            for i in range(start, start + len(pts) // 2 + 1):
                if sel(pts[i % len(pts)][2]): run.append('(%.3f,%.3f)' % pts[i % len(pts)][:2])
                else: break
            num.append(f'\\def\\ballArc{plane}{tag}{{{" ".join(run)}}}')
    summary['methane_ball'] = {'radius_cm': BALL_R, 'schematic': 'octahedron dual to the cube of third-turns, exact in Rodrigues coordinates',
                               'quarter_turn_true_radius': round(math.pi / 2 * bs, 4), 'drawn_radius': round(face_centre, 4)}
    # the loop: a third of a turn about the bond to hydrogen 1 returns X0 to its position with 2, 3, 4 cycled
    Rg = G.rot_axis(D[0], 2 * math.pi / 3)
    Xg = G.rotate(Rg, CH4)
    where = [min(range(4), key=lambda j: G.config_distance([Xg[i]], [CH4[j]])) for i in range(4)]   # nucleus i now sits where j was
    assert all(G.config_distance([Xg[i]], [CH4[where[i]]]) < 1e-9 for i in range(4)) and sorted(where) == [0, 1, 2, 3] and where[0] == 0
    summary['methane_loop'] = {'axis': 'C-H1', 'turn_deg': 120.0, 'nucleus_i_sits_where_j_was': [w + 1 for w in where]}
    # Td(M) = S4 on the four protons, classes (size, cycles, odd?, starred?): E, 3-cycles, double transpositions,
    # 4-cycles (starred), transpositions (starred); chi_spin = 2^cycles; chi_stat = sign; chi_pm = parity^star
    cls = [(1, 4, 1, 0), (8, 2, 1, 0), (3, 2, 1, 0), (6, 1, -1, 1), (6, 3, -1, 1)]
    TD = {'A1': [1, 1, 1, 1, 1], 'A2': [1, 1, 1, -1, -1], 'E': [2, -1, 2, 0, 0], 'T1': [3, 0, -1, 1, -1], 'T2': [3, 0, -1, -1, 1]}
    chi_spin = [2**c[1] for c in cls]
    spin_td = {k: int(round(sum(n*chi_spin[i]*v[i] for i, (n, _, _, _) in enumerate(cls))/24)) for k, v in TD.items()}
    weights = {}
    for par, tag in ((1, 'Even'), (-1, 'Odd')):
        total = [c[2]*(par if c[3] else 1) for c in cls]          # chi_stat chi_pm on each class
        weights[tag] = {k: int(round(sum(n*v[i]*chi_spin[i]*total[i] for i, (n, _, _, _) in enumerate(cls))/24)) for k, v in TD.items()}
    # the displacement representation on the fifteen Cartesian coordinates: chi(h) = (nuclei fixed by h) x tr(+-R_h^-1),
    # the sign from the star (the equivalent rotations: identity, third-turn, half-turn, quarter-turn, half-turn about
    # a cube edge; the starred classes act on displacements as improper operations)
    fixed = [5, 2, 1, 1, 3]
    trace = [3, 0, -1, -(1 + 0), -(1 - 2)]                      # tr(+-R): 1 + 2 cos theta, negated on the starred classes
    chi_3n = [f * t for f, t in zip(fixed, trace)]
    gamma_3n = {k: int(round(sum(n*c*v[i] for i, (n, _, _, _), c in zip(range(5), cls, chi_3n))/24)) for k, v in TD.items()}
    gamma_vib = dict(gamma_3n); gamma_vib['T1'] -= 1; gamma_vib['T2'] -= 1   # minus rotations (T1) and translations (T2)
    assert sum(m * TD[k][0] for k, m in gamma_3n.items()) == 15 and sum(m * TD[k][0] for k, m in gamma_vib.items()) == 9
    # D^J restricted to H through h -> R_h^-1 (the classes' turning angles), and the physical states per J by parity
    angles = [0.0, 120.0, 180.0, 90.0, 180.0]
    def chi_J(J, deg):
        if deg == 0.0: return 2 * J + 1
        t = math.radians(deg); return math.sin((2 * J + 1) * t / 2) / math.sin(t / 2)
    j_table = []
    for J in range(0, 7):
        ch = [chi_J(J, a) for a in angles]
        mult = {k: int(round(sum(n * c * v[i] for i, (n, _, _, _), c in zip(range(5), cls, ch)) / 24)) for k, v in TD.items()}
        assert sum(m * TD[k][0] for k, m in mult.items()) == 2 * J + 1
        j_table.append({'J': J, 'species': {k: m for k, m in mult.items() if m},
                        'even': sum(m * weights['Even'][k] for k, m in mult.items()), 'odd': sum(m * weights['Odd'][k] for k, m in mult.items())})
    # under the proper rotations T = A4: classes E, 3 C2, 4 C3, 4 C3' with chi_spin 16, 4, 4, 4; A, 1E, 2E, T
    w3 = cmath.exp(2j*math.pi/3)
    TT = {'A': [1, 1, 1, 1], '1E': [1, 1, w3, w3**2], '2E': [1, 1, w3**2, w3], 'T': [3, -1, 0, 0]}
    tcls = [1, 3, 4, 4]; tspin = [16, 4, 4, 4]
    spin_t = {k: int(round((sum(n*s*v.conjugate() for n, s, v in zip(tcls, tspin, map(complex, ch)))/12).real)) for k, ch in TT.items()}
    isomers = [('A', 'A', 1), ('1E', '2E', 1), ('2E', '1E', 1), ('T', 'T', 3)]      # (Gamma_rot, Gamma_nuc, d): Gamma_rot x Gamma_nuc contains A
    kernel_index = {k: int(round(12/sum(n for n, v in zip(tcls, ch) if abs(complex(v) - ch[0]) < 1e-9))) for k, ch in TT.items()}
    summary['methane_rigid'] = {'spin_Td': spin_td, 'weights': weights, 'spin_T': spin_t,
                                'isomers': [(r, nu, d, spin_t[nu]) for r, nu, d in isomers],
                                'entangled_fraction': [3*spin_t['T'], 16], 'monodromy_orders': kernel_index,
                                'gamma_3N': gamma_3n, 'gamma_vib': gamma_vib, 'J_table': j_table}

    # ---- methylamine
    X0 = G.methylamine(0.0, ETA0)
    summary['methylamine_frame'] = {k: round(v, 5) for k, v in G.MLA_FRAME.items()}
    summary['methylamine_X0'] = [[round(v, 4) for v in c] for c in X0]
    summary['methylamine_X0_Xm'] = [round(v, 12) for v in G.mass_moment(X0, G.MLA_MASSES)]
    emit_pic(mol, 'mol3d-mla-ref', pic_code_rods(X0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, labels=True))
    # figure 5 turns its glyphs in the page so that the projected half-turn axis of b (the body y axis) stands upright
    ax = G.project([(0.0, 1.0, 0.0)], CAM_MLA)[0]
    num.append(f'\\def\\mlaTurn{{{-math.degrees(math.atan2(-ax[0], ax[1])):.2f}}}')
    emit_pic(mol, 'mol3d-water-X', pic_code_rods(wX, G.WATER_ELEMENTS, bonds_w, CAM_FLAT, labels=True))
    RwX = G.rotate(Rw, wX)
    summary['water_matrices'] = {'X': [[round(v, 4) for v in c] for c in wX], 'RX': [[round(v, 4) for v in c] for c in RwX]}
    emit_pic(mol, 'mol3d-water-RX', pic_code_rods(G.rotate(Rw, wX), G.WATER_ELEMENTS, bonds_w, CAM_FLAT, labels=True))
    emit_pic(mol, 'mol3d-water-minusX', pic_code_rods(G.invert_config(wX), G.WATER_ELEMENTS, bonds_w, CAM_FLAT, labels=True))
    for k, (d, th) in enumerate(samples):
        emit_pic(mol, f'mol3d-water-sample-{k+1}', pic_code_rods(G.water(d, th), G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    # b X0 and R_b X0
    bX0 = G.apply_perm_inversion(X0, Q.b[0], Q.b[1])
    Rb = G.rot_y(math.pi)
    res_b = G.config_distance(bX0, G.rotate(Rb, X0))
    summary['b_equals_rotation_residual'] = res_b
    emit_pic(mol, 'mol3d-mla-bref', pic_code_rods(bX0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, labels=True))
    summary['mirror_plane_normal_is_Rb_axis'] = True
    # shapes along the t path at fixed orientation: tau = pi/3 (eclipsed) and 2 pi/3 (the version tH) (figure 6d)
    summary['tau_shapes_deg'] = [0, 60, 120]
    for deg in (0, 60, 120):
        emit_pic(mol, f'newman-tau-{deg}', pic_code_newman(G.methylamine(math.radians(deg), ETA0), G.MLA_ELEMENTS, twist_deg=(20.0 if deg == 60 else 0.0)))
    emit_pic(mol, 'newman-tuH', pic_code_newman(G.methylamine(-math.pi/3, -ETA0), G.MLA_ELEMENTS))
    # fragment loops in page coordinates (figure 4)
    P0 = G.project(X0, CAM_MLA)
    for name, members in (('methyl', [1, 2, 3, 6]), ('amino', [4, 5, 7])):
        pts = rounded_loop([(P0[i-1][0], P0[i-1][1]) for i in members], 0.78)
        num.append(f'\\def\\fragloop{name}{{' + ' '.join(f'({fmt(x)},{fmt(y)})' for x, y in pts) + '}')
    # version labels (figure 5c): position j carries label g(j)
    versions = version_positions()
    vlines = []
    for cname, gname, ta, et in versions:
        g = Q.E(Q.N)
        for ch in gname.replace('^2', '2'):
            pass
        # rebuild g from its name
        word = {'E': [], 't': ['t'], 't^2': ['t', 't'], 'u': ['u'], 'tu': ['t', 'u'], 't^2u': ['t', 't', 'u']}[gname]
        for w in word:
            g = Q.mul(g, {'t': Q.t, 'u': Q.u}[w])
        labels = [g[0][j] + 1 for j in range(7)]
        vlines.append(f'{{{cname}}}/{{{gname}}}/{"/".join(str(l) for l in labels)}')
    num.append('\\def\\versionList{' + ','.join(vlines) + '}')
    summary['versions'] = [(c, g, round(ta/math.pi, 6), round(et, 5)) for c, g, ta, et in versions]
    actions, worst = derive_family_actions()
    summary['family_actions'] = actions; summary['family_action_residual'] = worst
    # H-invariant cell U = [-pi/3, pi/3] x [0, eta0] and its five translates under the derived actions
    acts = {'t': (2/3, 1), 'u': (1.0, -1), 'tu': (-1/3, -1), 't^2': (4/3, 1), 't^2u': (1/3, -1)}
    cells = [('H', -1/3, 1/3, 0, 1)]
    for gname, (shift, sgn) in acts.items():
        lo, hi = -1/3 + shift, 1/3 + shift
        elo, ehi = (0, 1) if sgn > 0 else (-1, 0)
        # reduce to (-1, 1] in units of pi, splitting at the seam
        pieces = []
        lo2 = ((lo + 1) % 2) - 1; hi2 = lo2 + (hi - lo)
        if hi2 <= 1 + 1e-9:
            pieces.append((lo2, hi2))
        else:
            pieces.append((lo2, 1.0)); pieces.append((-1.0, hi2 - 2))
        for a, b in pieces:
            cells.append((gname + 'H', a, b, elo, ehi))
    num.append('\\def\\chartCells{' + ','.join(f'{{{n}}}/{a:.4f}/{b:.4f}/{e1}/{e2}' for n, a, b, e1, e2 in cells) + '}')
    summary['chart_cells'] = cells

    # normal frame and T_tu (figure 10)
    nf = normal_frame_data()
    summary['normal_frame'] = {'gram': nf['gram'], 'ortho': nf['ortho'], 'M': nf['M'],
                               'resid': nf['resid'], 'cfg_res': nf['cfg_res'],
                               'tangent_dim': nf['tangent_dim']}
    amp = 0.42        # actual normal displacement q (angstrom) used for the displaced drawing
    arrow_gain = 3.0  # arrows are drawn 2.5 x longer than the displacement; declared in the figure
    def arrows_for(e):
        Pe = [(G.dot(v, CAM_MLA[0]), G.dot(v, CAM_MLA[1])) for v in e]
        return [(arrow_gain*amp*px, arrow_gain*amp*py) if math.hypot(px, py)*amp > 0.04 else None for px, py in Pe]
    CAM_SIDE = G.camera(azimuth_deg=75.0, elevation_deg=18.0)   # the C-N axis lies in the page: both normal displacements at full length;
    # azimuth 75 (not 90) so that no methyl hydrogen sits on the line of sight through the carbon
    def arrows_side(e, longest=1.1):
        # side views: each direction's arrows scaled so that its largest arrow is `longest` angstrom on the page
        # (the direction is the content; a common gain cannot serve a hydrogen twist and a heavy-atom stretch); declared
        Pe = [(G.dot(v, CAM_SIDE[0]), G.dot(v, CAM_SIDE[1])) for v in e]
        big = max(math.hypot(px, py) for px, py in Pe)
        g = longest / big
        return [(g*px, g*py) if math.hypot(px, py) > 0.12*big else None for px, py in Pe]
    emit_pic(mol, 'mol3d-mla-ref-eminus-side', pic_code_rods(X0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_SIDE, arrows=arrows_side(nf['e_minus'])))
    emit_pic(mol, 'mol3d-mla-ref-eplus-side', pic_code_rods(X0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_SIDE, arrows=arrows_side(nf['e_plus'])))
    Xdisp = [G.add(x, G.scale(v, amp)) for x, v in zip(X0, nf['e_minus'])]
    Rdisp = G.rot_axis((0.3, 1.0, 0.25), math.radians(55.0))
    emit_pic(mol, 'mla-disp', pic_code_rods(Xdisp, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    # the undisplaced ghost beneath it (figure 10): only the atoms that visibly move and the rods that touch them; an
    # atom that moves less than a quarter of a hydrogen radius on the page would show only as a sliver beside its copy
    P0, Pd = G.project(X0, CAM_MLA), G.project(Xdisp, CAM_MLA)
    shift = [math.hypot(Pd[k][0] - P0[k][0], Pd[k][1] - P0[k][1]) for k in range(len(X0))]
    moved = {k for k in range(len(X0)) if shift[k] > 0.05}
    assert moved == {3, 4} and min(shift[k] for k in moved) > 0.2 and max(shift[k] for k in range(len(X0)) if k not in moved) < 0.04, shift
    emit_pic(mol, 'mla-disp-ghost', pic_code_rods(X0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, only=moved))
    summary['mla_disp_page_shift'] = [round(v, 4) for v in shift]
    emit_pic(mol, 'mla-rot-disp', pic_code_rods(G.rotate(Rdisp, Xdisp), G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))

    # generic X and bX (figure 11b): X = X0(0.22, eta0) + small normal displacement
    Xg = G.methylamine(0.22, ETA0)
    eg, _ = G.normal_vector([(0,0,0)]*3 + [(0,0,1.0), (0,0,-1.0)] + [(0,0,0)]*2, 0.22, ETA0)
    Xg = [G.add(x, G.scale(v, 0.30)) for x, v in zip(Xg, eg)]
    bXg = G.apply_perm_inversion(Xg, Q.b[0], Q.b[1])
    R, res = G.kabsch_rotation(Xg, bXg, G.MLA_MASSES)
    summary['generic_X_bX_best_rotation_residual'] = res
    summary['generic_X_Xm'] = [round(v, 12) for v in G.mass_moment(Xg, G.MLA_MASSES)]
    emit_pic(mol, 'mol3d-mla-X', pic_code_rods(Xg, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, labels=True))
    emit_pic(mol, 'mol3d-mla-bX', pic_code_rods(bXg, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, labels=True))
    # free-action check near the reference: min_g ||gX0 - X0|| over g != E in G12
    dmin = min(G.config_distance(G.apply_perm_inversion(X0, g[0], g[1]), X0) for g in Q.G12 if g != Q.E(Q.N))
    summary['min_displacement_by_nontrivial_G12_element'] = dmin

    # continuation results
    tb = Q.mul(Q.t, Q.b); bt = Q.mul(Q.b, Q.t); t2b = Q.mul(Q.mul(Q.t, Q.t), Q.b)
    summary['bt_equals_t2b'] = (bt == t2b); summary['tb_ne_bt'] = (tb != bt)
    summary['btb_equals_tinv'] = (Q.mul(Q.b, Q.mul(Q.t, Q.b)) == Q.inv(Q.t))

    # ---- KRb (figure 4b)
    K = G.krb_schematic()
    Sk, Gin, GK, GRb, p12, p34, Es = Q.krb_groups()
    summary['krb_orders'] = [len(Sk), len(Gin), len(GK), len(GRb)]
    summary['krb_any_two_generate'] = all(len(Q.close(list(a | b), 4)) == 8 for a, b in [(Gin, GK), (Gin, GRb), (GK, GRb)])
    # the three channel groups meet exactly in <E*> and cover S (figure 4c, nested outlines)
    core = Gin & GK & GRb
    assert core == (Gin & GK) == (Gin & GRb) == (GK & GRb) == {Q.E(4), Es}, 'channel groups must meet in <E*>'
    assert Gin | GK | GRb == Sk, 'channel groups must cover S'
    def krb_name(g):
        perm, star = g
        cyc = []
        if perm[0] == 1: cyc.append('(12)')
        if perm[2] == 3: cyc.append('(34)')
        base = ''.join(cyc) if cyc else 'E'
        return base + ('^*' if star else '')
    def elist(S_):
        return ', '.join(krb_name(g) for g in sorted(S_, key=lambda g: (bool(g[1]), g[0])))
    summary['krb_core'] = elist(core); summary['krb_lobes'] = {'in': elist(Gin - core), 'K': elist(GK - core), 'Rb': elist(GRb - core)}
    for name, lines_ in (('in', [(1, 3), (2, 4)]), ('complex', [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]), ('out', [(1, 2), (3, 4)])):
        emit_pic(mol, f'mol3d-krb-{name}', pic_code_grouped(K, G.KRB_ELEMENTS, lines_, CAM_FLAT))
    # every subgroup of S (order 8): generated by at most three elements
    allsubs = set()
    els = sorted(Sk)
    for a in els:
        for b in els:
            for c in els:
                allsubs.add(frozenset(Q.close([a, b, c], 4)))
    allsubs = sorted(allsubs, key=lambda K_: (len(K_), sorted(K_)))
    assert len(allsubs) == 16, len(allsubs)
    named = {frozenset([Q.E(4)]): 'E', frozenset(core): '\\langle E^*\\rangle', frozenset(Gin): 'G_{\\mathrm{in}}',
             frozenset(GK): 'G_{\\mathrm K}', frozenset(GRb): 'G_{\\mathrm{Rb}}', frozenset(Sk): 'S'}
    rows = {}
    for K_ in allsubs: rows.setdefault(len(K_), []).append(K_)
    # centre the named subgroups in their rows
    for L, row in rows.items():
        row.sort(key=lambda K_: (0 if K_ in named else 1, sorted(K_)))
        n_ = len(row); named_in = [K_ for K_ in row if K_ in named]; others = [K_ for K_ in row if K_ not in named]
        order_ = []
        mid = (n_ - 1)/2
        slots = sorted(range(n_), key=lambda i: abs(i - mid))
        placed = {}
        for K_, slot in zip(named_in + others, slots): placed[slot] = K_
        rows[L] = [placed[i] for i in range(n_)]
    key_ = {}; nodes = []
    ylev = {1: 0.0, 2: 1.0, 4: 2.0, 8: 3.0}
    for L, row in rows.items():
        for i, K_ in enumerate(row):
            key_[K_] = len(nodes)
            x = (i + 0.5)/len(row)
            nodes.append((key_[K_], x, ylev[L], 1 if K_ in named else 0, named.get(K_, '')))
    cov_ = []
    for A in allsubs:
        for B in allsubs:
            if A < B and len(B) == 2*len(A):
                cov_.append((key_[A], key_[B], 1 if (A in named and B in named) else 0))
    num.append('\\def\\krbNodes{' + ','.join(f'{i}/{x:.4f}/{y:.1f}/{m}/{{{nm}}}' for i, x, y, m, nm in nodes) + '}')
    idx = {named[K_]: key_[K_] for K_ in named}
    core_name = '\\langle E^*\\rangle'; in_name = 'G_{\\mathrm{in}}'; k_name = 'G_{\\mathrm K}'; rb_name = 'G_{\\mathrm{Rb}}'
    num.append('\\def\\krbIdxE{%d}\\def\\krbIdxCore{%d}\\def\\krbIdxIn{%d}\\def\\krbIdxK{%d}\\def\\krbIdxRb{%d}\\def\\krbIdxS{%d}'
               % (idx['E'], idx[core_name], idx[in_name], idx[k_name], idx[rb_name], idx['S']))
    num.append('\\def\\krbCovers{' + ','.join(f'{a}/{b}/{m}' for a, b, m in cov_) + '}')
    summary['krb_subgroups'] = len(allsubs); summary['krb_covers'] = len(cov_)

    # ---- subgroup interval (figure 5b)
    subs, names, cov = Q.interval_data()
    order_of = {K: len(K)//2 for K in subs}     # |K/H|
    # drawing x-slots per node (design choice); y = log2 of the version count |K/H|
    xslot = {'H': 0.0, '\\langle H,u\\rangle': 4.6, '\\langle H,E^*\\rangle': -0.4, '\\langle H,u^*\\rangle': -2.9,
             'G_6': 2.6, '\\langle H,u,E^*\\rangle': -2.6, 'G_{12}': 3.2, '\\langle G_6,E^*\\rangle': -0.6,
             '\\langle G_6,u^*\\rangle': -4.4, 'B': 0.0}
    ordered = sorted(subs, key=lambda K: (len(K), names[K]))
    key = {names[K]: i for i, K in enumerate(ordered)}
    num.append('\\def\\hasseNodes{' + ','.join(
        f'{key[names[K]]}/{xslot[names[K]]}/{math.log2(len(K)//2):.4f}/{len(K)//2}/{{{names[K]}}}' for K in ordered) + '}')
    # the generator added along each cover of [H, B], and a minimal generator form of every subgroup (b always first)
    from itertools import combinations
    cands = [('t', Q.t), ('u', Q.u), ('E^*', Q.Estar), ('u^*', Q.mul(Q.u, Q.Estar)), ('t^*', Q.mul(Q.t, Q.Estar))]
    labels = []
    for a_, b_ in cov:
        lab = '?'
        for nm, g in cands:
            if frozenset(Q.close(list(a_) + [g], Q.N)) == frozenset(b_):
                lab = nm; break
        assert lab != '?', 'cover without a single added generator'
        labels.append((key[names[a_]], key[names[b_]], lab))
    num.append('\\def\\hasseCoverLabels{' + ','.join(f'{a}/{b}/{{{nm}}}' for a, b, nm in labels) + '}')
    gens = {}
    chain_forms = {'H': [], 'G_6': ['t'], 'G_{12}': ['t', 'u'], 'B': ['t', 'u', 'E^*']}
    cand_by_name = dict(cands)
    for K_ in subs:
        found = None
        if names[K_] in chain_forms:
            form = chain_forms[names[K_]]
            assert frozenset(Q.close([Q.b] + [cand_by_name[nm] for nm in form], Q.N)) == frozenset(K_), names[K_]
            found = ','.join(['b'] + form)
        else:
            for r in range(0, 4):
                for combo in combinations(cands, r):
                    if frozenset(Q.close([Q.b] + [g for _, g in combo], Q.N)) == frozenset(K_):
                        found = ','.join(['b'] + [nm for nm, _ in combo]); break
                if found: break
        assert found, 'no generator form'
        gens[key[names[K_]]] = found
    summary['bond_interval_generators'] = {names[K_]: gens[key[names[K_]]] for K_ in subs}
    chainnames = {'H': 'H', 'G_6': 'G_6', 'G_{12}': 'G_{12}', 'B': 'B'}
    labs = []
    for K_ in subs:
        k = key[names[K_]]; gform = '\\langle ' + gens[k] + '\\rangle'
        labs.append((k, (chainnames[names[K_]] + '=' + gform) if names[K_] in chainnames else gform))
    num.append('\\def\\hasseLabels{' + ','.join(f'{k}/{{{l}}}' for k, l in sorted(labs)) + '}')
    # the whole interval [H, S]: every subgroup containing H, with covers; the bond interval is a sub-poset
    big = Q.interval(Q.H, Q.S, Q.N)
    bigcov = Q.covers(big)
    bigsorted = sorted(big, key=lambda K: (len(K), sorted(K)))
    bigkey = {K: i for i, K in enumerate(bigsorted)}
    inbond = {bigkey[K] for K in big if K <= Q.B}
    levels = sorted(set(len(K) for K in big))
    lev = {K: levels.index(len(K)) for K in big}
    yof = {K: math.log(len(K)/2)/math.log(120) for K in big}      # height by log of the version count
    below = {K: [a for a, b in bigcov if b == K] for K in big}
    above = {K: [b for a, b in bigcov if a == K] for K in big}
    rows = {L: [K for K in bigsorted if len(K) == L] for L in levels}

    def layout(pull, median, sweeps):
        """Layered placement: y by level; x by repeated barycentre (or median) ordering of each level from its
        neighbours, alternating the sweep direction, nodes spread evenly in each level; the bond interval pulled
        left by `pull` so that [H, B] reads as one side of the picture."""
        xpos = {K: 0.5 for K in big}
        for sweep in range(sweeps):
            for L in (levels if sweep % 2 == 0 else levels[::-1]):
                row = rows[L]
                def key(K):
                    nb = (below[K] if sweep % 2 == 0 else above[K]) or (below[K] + above[K])
                    vals = sorted(xpos[a] for a in nb)
                    if not vals:
                        v = xpos[K]
                    elif median:
                        v = vals[len(vals)//2] if len(vals) % 2 else (vals[len(vals)//2 - 1] + vals[len(vals)//2])/2
                    else:
                        v = sum(vals)/len(vals)
                    return v + (-pull if K <= Q.B else 0.0)
                order = sorted(row, key=lambda K: (key(K), bigkey[K]))
                for i, K in enumerate(order):
                    xpos[K] = (i + 0.5)/len(order)
        return xpos

    def crossings(xpos):
        """Number of pairs of cover segments that cross in the drawing (shared endpoints do not count)."""
        def seg(a, b): return (xpos[a], yof[a], xpos[b], yof[b])
        def cross(p, q):
            (x1, y1, x2, y2), (x3, y3, x4, y4) = p, q
            def orient(ax, ay, bx, by, cx, cy): return (bx - ax)*(cy - ay) - (by - ay)*(cx - ax)
            o1, o2 = orient(x1, y1, x2, y2, x3, y3), orient(x1, y1, x2, y2, x4, y4)
            o3, o4 = orient(x3, y3, x4, y4, x1, y1), orient(x3, y3, x4, y4, x2, y2)
            return o1*o2 < 0 and o3*o4 < 0
        n = 0
        for i, (a1, b1) in enumerate(bigcov):
            for a2, b2 in bigcov[i+1:]:
                if len({a1, b1, a2, b2}) < 4: continue
                if cross(seg(a1, b1), seg(a2, b2)): n += 1
        return n

    best = None
    for pull in (0.0, 0.15, 0.25, 0.35):
        for median in (False, True):
            for sweeps in (8, 16, 24):
                xp = layout(pull, median, sweeps)
                score = (crossings(xp), -pull)      # fewest crossings; among equals, the bond interval further left
                if best is None or score < best[0]:
                    best = (score, xp, (pull, median, sweeps))
    xpos = best[1]
    summary['interval_HS_layout'] = {'crossings': best[0][0], 'pull': best[2][0], 'median': best[2][1], 'sweeps': best[2][2]}
    num.append('\\def\\bigNodes{' + ','.join(f'{bigkey[K]}/{xpos[K]:.4f}/{yof[K]:.4f}/{1 if bigkey[K] in inbond else 0}/{{{names.get(K, "")}}}' for K in bigsorted) + '}')
    num.append('\\def\\bigCovers{' + ','.join(f'{bigkey[a]}/{bigkey[b]}/{1 if (a <= Q.B and b <= Q.B) else 0}' for a, b in bigcov) + '}')
    chain_keys = {names[K]: bigkey[K] for K in big if K in names and names[K] in ('H', 'G_6', 'G_{12}', 'B')}
    num.append('\\def\\bigChain{' + ','.join(str(chain_keys[n]) for n in ('H', 'G_6', 'G_{12}', 'B')) + '}')
    chain_seq = [chain_keys[n] for n in ('H', 'G_6', 'G_{12}', 'B')]
    num.append('\\def\\bigChainCovers{' + ','.join(f'{a}/{b}' for a, b in zip(chain_seq, chain_seq[1:])) + '}')
    summary['interval_HS'] = {'subgroups': len(big), 'covers': len(bigcov), 'in_bond_interval': len(inbond),
                              'orders': sorted(len(K) for K in big)}
    summary['interval'] = {'subgroups': len(subs), 'covers': len(cov),
                           'versions': sorted(order_of.values())}

    # ---- species and spin weights (figures 7, 8)
    cnames, act_t, act_b = Q.coset_action_table()
    summary['coset_action'] = {'t': act_t, 'b': act_b}
    w = Q.spin_weights()
    summary['spin_weights'] = {str(k): v for k, v in w.items()}
    num.append(f'\\def\\spinAone{{{w[1]["A1"]}}}\\def\\spinAtwo{{{w[1]["A2"]}}}\\def\\spinE{{{w[1]["E"]}}}\\def\\spinDim{{{w[1]["A1"]+w[1]["A2"]+2*w[1]["E"]}}}')   # figure 8, the proton factor
    # multiplicities of species in C[G6/H]: permutation character (3, 0, 1)
    perm_char = {'E': 3, 't': 0, 'b': 1}
    local = {n: int(round(sum(Q.G6_CHARS[n][c]*perm_char[c]*{'E': 1, 't': 2, 'b': 3}[c] for c in 'Etb')/6)) for n in Q.G6_CHARS}
    summary['local_copies'] = local

    # ---- Gaussian shares (figure 9)
    with open(os.path.join(DATA, 'gaussian-shares.dat'), 'w') as f:
        f.write('x wA1 wE wA2\n')
        f.write('0.00 0.333333 0.666667 0\n')   # the limit: no overlap, shares 1/3 and 2/3
        for i in range(2, 201):
            x = i/100
            c = math.exp(-1/(8*x*x))
            f.write(f'{x:.2f} {(1+2*c)/3:.6f} {2*(1-c)/3:.6f} 0\n')
    # figure 9's three rows: Delta/d, the packet width in degrees for d = 120 degrees, the A1 share
    rows = [(r, (1 + 2*math.exp(-1/(8*float(r)**2)))/3) for r in (Fraction(1, 6), Fraction(1, 3), Fraction(2, 3))]
    num.append('\\def\\packetRows{' + ','.join(f'{{{r}}}/{float(r)*120:g}/{w:.4f}' for r, w in rows) + '}')
    num.append('\\def\\packetMarks{' + ' '.join(f'({float(r):.4f},{w:.4f})' for r, w in rows) + '}')
    x_show = 1/3; c_show = math.exp(-1/(8*x_show*x_show))
    summary['gaussian_shown'] = {'Delta_over_d': x_show, 'c': c_show, 'wA1': (1+2*c_show)/3, 'wE': 2*(1-c_show)/3}

    # ---- component functions (figure 2): a physical choice on the delta = 0 slice.  There X P_sigma = R X with R the
    # in-plane half-turn, so a rotation-invariant state obeys Psi(X) = chi_stat(sigma) sigma.Psi(X) = -(b xi + a eta),
    # i.e. b = -a: the proton singlet (xi - eta) times a scalar f(theta).  f is an illustrative Gaussian in theta.
    COMP_CENTRE, COMP_WIDTH = 104.5, 14.0   # degrees
    num.append(f'\\def\\compCentre{{{COMP_CENTRE}}}\\def\\compWidth{{{COMP_WIDTH}}}')
    def comp_a(th): return math.exp(-((th-COMP_CENTRE)/COMP_WIDTH)**2/2)
    def comp_b(th): return -comp_a(th)
    with open(os.path.join(DATA, 'component-functions.dat'), 'w') as f:
        f.write('theta a b n2\n')
        for i in range(0, 101):
            th = 80 + i
            f.write(f'{th} {comp_a(th):.5f} {comp_b(th):.5f} {comp_a(th)**2+comp_b(th)**2:.5f}\n')
    for k, (d, th) in enumerate(samples):
        tag = 'one two three'.split()[k]
        num.append(f'\\def\\sampleTheta{tag}{{{th:g}}}\\def\\sampleA{tag}{{{comp_a(th):.2f}}}\\def\\sampleB{tag}{{{comp_b(th):+.2f}}}')
    summary['sample_components'] = [(th, round(comp_a(th), 4), round(comp_b(th), 4)) for d, th in samples]

    # ---- momentum labels (figure 12)
    rho = 0.5          # the symmetric frame: each turn of the frame is half the torsion (Mellor, Yurchenko, Mant, Jensen 2019)
    num.append(f'\\def\\rhoIll{{{rho}}}')
    rows = []
    for Kq in (-2, -1, 0, 1, 2):
        rows.append(f'{Kq}/' + '{' + ','.join(f'{m + rho*Kq:.1f}' for m in range(-3, 4)) + '}')
    num.append('\\def\\kappaRows{' + ','.join(rows) + '}')
    tms = torsion_mode_species()
    num.append('\\def\\torsionSpecies{' + ','.join(f'{m}/{c}/{s if s else "none"}' for m, c, s in tms) + '}')
    summary['torsion_species'] = tms
    summary['kappa_rule'] = 'kappa = m + rho K from (r,tau+2pi,q) ~ (r R_z(-2 pi rho),tau,q) and D^J_MK(r R_z(omega)) = e^{-iK omega} D^J_MK(r)'

    with open(os.path.join(DATA, 'molecules.tex'), 'w') as f:
        f.write('\n'.join(mol) + '\n')
    with open(os.path.join(DATA, 'numbers.tex'), 'w') as f:
        f.write('\n'.join(num) + '\n')
    with open(os.path.join(DATA, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=1, default=float)
    print('wrote', DATA)
    print(json.dumps({k: summary[k] for k in ('b_equals_rotation_residual', 'family_action_residual',
          'normal_frame', 'generic_X_bX_best_rotation_residual', 'min_displacement_by_nontrivial_G12_element',
          'versions', 'spin_weights', 'local_copies', 'torsion_species', 'interval', 'krb_any_two_generate')}, indent=1, default=float))

if __name__ == '__main__':
    main()
