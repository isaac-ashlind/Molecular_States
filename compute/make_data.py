"""Generate the figure data in figures/data/ from compute/.

Run from the repository root:  python3 compute/make_data.py
numbers.tex and molecules.tex hold what the plates draw, summary.json what the
checks compare.  Everything is computed from the masses, the model parameters
and the group definitions, with characters and counts in exact arithmetic;
the character tables of Td(M) and T and the drawing choices are typed in and
marked as such.  checks/ recompute what is stated and drawn.
"""
import itertools, json, math, os, sys
from fractions import Fraction
sys.path.insert(0, os.path.dirname(__file__))
import geometry as G
import groups as Q

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'figures', 'data')
os.makedirs(DATA, exist_ok=True)

IOTA0 = G.MLA_FRAME['iota0']
CAM_MLA = G.camera(azimuth_deg=50.0, elevation_deg=25.0)
CAM_FLAT = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))   # page = xy-plane

def fmt(x):
    s = f'{x:.4f}'
    return '0.0000' if s == '-0.0000' else s

# --------------------------------------------------------------- pic emit ---

def pic_code_rods(X, elements, bonds, cam, arrows=None, labels=False, only=None):
    """TikZ pic code for the molecule primitives (\\molatom, \\molrodj, \\molarrow, \\mollabel) with correct layering.

    Items (atoms and rods) are painted back to front.  A rod is painted just after its farther atom
    and before its nearer atom.  It starts where the bond cylinder leaves the far sphere (the junction
    circle, projected as a half-ellipse computed in TeX from the atom radius) and vanishes under the
    silhouette of the near sphere.  Rod depth = depth of the far atom + a small epsilon.
    `only`, a set of atom indices, restricts the drawing to those atoms and the rods that touch them.
    """
    right, up, out = cam
    P = G.project(X, cam)
    n = len(X)
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
            assert dn > 1e-6, 'an end-on bond'
            rad = RAD_MACRO.get(elements[far], '\\radX')
            lines.append(f'\\molrodj{{{fmt(P[far][0])}}}{{{fmt(P[far][1])}}}{{{fmt(P[near][0])}}}{{{fmt(P[near][1])}}}'
                         f'{{{fmt(ux/dn)}}}{{{fmt(uy/dn)}}}{{{fmt(uz)}}}{{{rad}}}')
    if arrows is not None:
        rad_a = {'H': 0.20, 'C': 0.36, 'N': 0.36, 'O': 0.36, 'K': 0.40, 'Rb': 0.46}   # \\radH, \\radX, ... of primitives.tex
        for k in range(n):
            if arrows[k] is not None:
                dx, dy = arrows[k]
                # the shaft starts on the atom's silhouette, not at its center (visible on a light disc)
                L = math.hypot(dx, dy); r = rad_a[elements[k]]
                assert L > r + 0.15, 'an arrow shorter than its atom'
                ux, uy = dx/L, dy/L
                lines.append(f'\\molarrow{{{fmt(P[k][0] + r*ux)}}}{{{fmt(P[k][1] + r*uy)}}}{{{fmt(dx - r*ux)}}}{{{fmt(dy - r*uy)}}}')
    if labels:
        # label direction: away from the mean page position of the bonded neighbors (or straight up)
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

def pic_code_newman(X, methyl=(1, 2, 3), amino=(4, 5), twist_deg=0.0):
    """Newman projection down the C-N axis (viewer on the N side): back carbon as a circle with its
    hydrogens on fine bonds from the rim; front nitrogen as a disc at the center with its hydrogens on
    rod bonds from the center.  Page coordinates are the molecular (x, y) in angstrom; the digit is the
    column label (version labels are supplied by the \\molA.. macros).  twist_deg turns the back set on
    the page by a small angle, the drawing convention for an eclipsed projection, not geometry."""
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
    lines.append('\\nmcenter')
    return lines

def exact_mean(total, order):
    """A character inner product: the integer total/order, asserted exact."""
    q, r = divmod(total, order)
    assert r == 0
    return q

def emit_pic(out, name, lines):
    out.append(f'\\tikzset{{pics/{name}/.style={{code={{%')
    for l in lines:
        out.append('  ' + l + '%')
    out.append('}}}')

def hull(points):
    pts = sorted(set(points))
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

    Each hull vertex is replaced by three points on the circle of radius
    margin about it (along the two edge normals and the bisector), so that
    plot[smooth cycle] through them hugs the hull with rounded corners.
    """
    h = hull(points)
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

def version_positions():
    """(coset, word, g, tau_g, iota_g) for the six coset representatives of G12/H, found among the 12 positions."""
    reps = [('H', [], 'E'), ('tH', ['t'], 't'), ('t^2H', ['t', 't'], 't^2'),
            ('uH', ['u'], 'u'), ('tuH', ['t', 'u'], 'tu'), ('t^2uH', ['t', 't', 'u'], 't^2u')]
    X0 = G.methylamine(0.0, IOTA0)
    out = []
    for cname, word, gname in reps:
        g = Q.E(Q.N)
        for w in word:
            g = Q.mul(g, {'t': Q.t, 'u': Q.u}[w])
        gX = G.apply_perm_inversion(X0, g[0], g[1])
        res, ta, io = min((G.kabsch_rotation(G.methylamine(ta, io), gX, G.MLA_MASSES)[1], ta, io)
                          for ta in [k*math.pi/3 for k in range(6)] for io in (IOTA0, -IOTA0))
        assert res < 1e-9, (cname, res)
        out.append((cname, gname, g, ta, io))
    return out

def normal_frame_data():
    """Two symmetry-adapted normal vectors at a = (0, iota0), the C-N stretch e_+ and the amino twist e_-, and their
    transformation under tu."""
    tau, iota = 0.0, IOTA0
    tau2, iota2 = tau - math.pi/3, -iota
    zero = (0.0, 0.0, 0.0)
    stretch = [zero]*7; stretch[5] = (0, 0, -1.0); stretch[6] = (0, 0, 1.0)
    twist = [zero]*7; twist[3] = (0, 0, 1.0); twist[4] = (0, 0, -1.0)
    e_plus, tangent = G.normal_vector(stretch, tau, iota)
    e_minus, _ = G.normal_vector(twist, tau, iota)
    ep_plus, _ = G.normal_vector(stretch, tau2, iota2)
    ep_minus, _ = G.normal_vector(twist, tau2, iota2)
    m = G.MLA_MASSES
    ortho = max(abs(G.mass_inner(e, tv, m)) for e in (e_plus, e_minus) for tv in tangent)
    # transformation under tu: e_a(a) P_tu = A sum_b e_b(a') M_ba, with A = R_z(pi)
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
    return dict(e_plus=e_plus, e_minus=e_minus, ortho=ortho, M=M, resid=resid, tangent_dim=len(tangent))

def torsion_mode_species():
    """Species of cos(m tau), sin(m tau) under G6 with t: tau -> tau + 2pi/3, b: tau -> -tau."""
    out = []
    for m in range(0, 7):
        if m == 0:
            out.append((0, 'A1', None)); continue
        # rep on (cos, sin): U_t = rotation by -2 pi m/3 (trace 2 cos(2 pi m/3) = 2 or -1, exactly), U_b = diag(1,-1)
        chi_t = 2 if m % 3 == 0 else -1; chi_b = 0; chi_E = 2
        mult = {name: exact_mean(chi['E']*chi_E + 2*chi['t']*chi_t + 3*chi['b']*chi_b, 6) for name, chi in Q.G6_CHARS.items()}
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
    bonds_w = [(1, 3), (2, 3)]
    Rw = G.rot_z(math.radians(40.0))
    # the shape grid of figure 1: delta ratios and angles, drawn mass-centered
    grid_d = [-0.5, -0.25, 0.0, 0.25, 0.5]; grid_t = [60.0, 90.0, 120.0, 150.0, 180.0]
    for i, d in enumerate(grid_d):
        for j, th in enumerate(grid_t):
            Xg = G.water(d, th)
            emit_pic(mol, f'mol3d-water-grid-{i}{j}', pic_code_rods(Xg, G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    # three sample configurations for figure 2
    samples = [(0.0, 90.0), (0.0, 110.0), (0.0, 140.0)]    # on the delta = 0 slice, in order of theta

    # ---- methane (the closing pages): the rigid tetrahedron in its C2 frame, the orientation ball SO(3) with the twelve
    # symmetry rotations and one cell, the vibrational species, the J ladder; the numbers against Albert et al.
    CH4 = G.methane_c2()
    # one camera for both panels, so that bond 1 on the molecule points along the arrow in the ball: G.camera(108, 21)
    # with the x and z components of its three vectors exchanged. G.camera builds a left-handed frame (a mirror image);
    # the exchange x <-> z is a mirror symmetry of the cube, the dual octahedron and the bond tetrahedron (it fixes bond
    # 1 and exchanges hydrogens 2 and 4), so the frame becomes right-handed and the molecule a true view with x to the
    # right, y up and z toward the viewer.
    def _swap_xz(v): return (v[2], v[1], v[0])
    CAM_BALL = tuple(_swap_xz(v) for v in G.camera(azimuth_deg=108.0, elevation_deg=21.0))
    assert G.dot(G.cross(CAM_BALL[0], CAM_BALL[1]), CAM_BALL[2]) > 0.999          # right-handed
    assert _swap_xz(tuple(CH4[0])) == tuple(CH4[0]) and _swap_xz(tuple(CH4[1])) == tuple(CH4[3])   # x <-> z fixes H1 and exchanges H2 and H4
    emit_pic(mol, 'mol3d-ch4-c2', pic_code_rods(CH4, G.CH4_ELEMENTS, G.CH4_BONDS, CAM_BALL, labels=True))
    front, back = [], []
    for name, v in (('X', (1.0, 0.0, 0.0)), ('Y', (0.0, 1.0, 0.0)), ('Z', (0.0, 0.0, 1.0))):
        px, py, pz = G.project([v], CAM_BALL)[0]
        num.append(f'\\def\\chAxis{name}x{{{px:.4f}}}\\def\\chAxis{name}y{{{py:.4f}}}')
        (front if pz >= 0 else back).append(name + '/1'); (back if pz >= 0 else front).append(name + '/-1')   # each half-axis: toward the viewer over the molecule, away from it behind
    num.append('\\def\\chAxisFront{' + ','.join(front) + '}\\def\\chAxisBack{' + ','.join(back) + '}')
    # the orientation space as the axis-angle ball (direction the axis, distance the angle, radius pi; the skin is the
    # half-turns, each point its own antipode). The twelve rotations of T: the center, the eight third-turns at 2pi/3
    # along the body diagonals, the three half-turns on the skin along the C2 axes. The cell of the center is drawn as
    # the schematic exact in Rodrigues coordinates: the octahedron dual to the cube spanned by the eight third-turns,
    # its vertices (the quarter-turns) on the cube's face centers, every edge straight; in the ball the true edges bow
    # outward (the caption says so). The cube itself is a reading aid: no edge between two points is intrinsic.
    BALL_R = 2.0                                                # cm for the angle pi: the ball's radius on the page
    bs = BALL_R / math.pi                                       # cm per radian
    D = [G.unit(CH4[k]) for k in range(4)]                      # bond directions = body diagonals
    outv = CAM_BALL[2]
    def emit_point(name, v, scale=1.0):
        px, py, _ = G.project([G.scale(v, scale)], CAM_BALL)[0]
        num.append(f'\\def\\{name}x{{{px:.4f}}}\\def\\{name}y{{{py:.4f}}}')
    num.append(f'\\def\\ballR{{{BALL_R:.4f}}}')
    for k, tag in enumerate('abcd'):
        emit_point('ballC' + tag, D[k], 2 * math.pi / 3 * bs)   # third-turn about bond k, +2pi/3: a cube vertex
        emit_point('ballD' + tag, D[k], -2 * math.pi / 3 * bs)  # the opposite turn
    axes = {'x': (1.0, 0.0, 0.0), 'y': (0.0, 1.0, 0.0), 'z': (0.0, 0.0, 1.0)}
    face_center = 2 * math.pi / 3 * bs / math.sqrt(3.0)          # the cube's half-edge: where its face centers sit on the axes
    for a_, v in axes.items():
        for sgn, t in ((1.0, 'p'), (-1.0, 'm')):
            if sgn * G.dot(v, outv) < 0: emit_point(f'ballS{a_}{t}', v, sgn * math.pi * bs)   # a half-turn, at its far side on the skin
            emit_point(f'ballV{a_}{t}', v, sgn * face_center)           # quarter-turn, drawn on the cube's face center: a vertex of the dual octahedron
    emit_point('ballFa', D[0], 2 * math.pi / 9 * bs)             # where the lift leaves the octahedron: the centroid of its face toward bond 1 (the dual octahedron's face plane x+y+z = c lies at c/sqrt3 on the diagonal)
    num.append('\\def\\ballSkinBack{' + ','.join(f'{a_}{t}' for a_, v in axes.items() for sgn, t in ((1.0, 'p'), (-1.0, 'm')) if sgn * G.dot(v, outv) < 0) + '}')   # each half-turn once, at its far side, drawn pale
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
    # under the front lines, like the far half-turns
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
    # the loop: a third of a turn about the bond to hydrogen 1 returns X0 to its position with 2, 3, 4 cycled
    Rg = G.rot_axis(D[0], 2 * math.pi / 3)
    Xg = G.rotate(Rg, CH4)
    where = [min(range(4), key=lambda j: G.config_distance([Xg[i]], [CH4[j]])) for i in range(4)]   # nucleus i now sits where j was
    assert all(G.config_distance([Xg[i]], [CH4[where[i]]]) < 1e-9 for i in range(4)) and sorted(where) == [0, 1, 2, 3] and where[0] == 0
    summary['methane_loop'] = {'axis': 'C-H1', 'nucleus_i_sits_where_j_was': [w + 1 for w in where]}
    # Td(M) = S4 on the four protons, classes (size, cycles, odd?, starred?): E, 3-cycles, double transpositions,
    # 4-cycles (starred), transpositions (starred); chi_spin = 2^cycles; chi_stat = sign; chi_pm = parity^star. The
    # character table is typed; checks/verify_methane.py and verify.g derive these numbers independently
    cls = [(1, 4, 1, 0), (8, 2, 1, 0), (3, 2, 1, 0), (6, 1, -1, 1), (6, 3, -1, 1)]
    TD = {'A1': [1, 1, 1, 1, 1], 'A2': [1, 1, 1, -1, -1], 'E': [2, -1, 2, 0, 0], 'T1': [3, 0, -1, 1, -1], 'T2': [3, 0, -1, -1, 1]}
    chi_spin = [2**c[1] for c in cls]
    spin_td = {k: exact_mean(sum(n*chi_spin[i]*v[i] for i, (n, _, _, _) in enumerate(cls)), 24) for k, v in TD.items()}
    weights = {}
    for par, tag in ((1, 'Even'), (-1, 'Odd')):
        total = [c[2]*(par if c[3] else 1) for c in cls]          # chi_stat chi_pm on each class
        weights[tag] = {k: exact_mean(sum(n*v[i]*chi_spin[i]*total[i] for i, (n, _, _, _) in enumerate(cls)), 24) for k, v in TD.items()}
    # the displacement representation on the fifteen Cartesian coordinates: chi(h) = (nuclei fixed by h) x tr(+-R_h^-1),
    # the sign from the star (the equivalent rotations: identity, third-turn, half-turn, quarter-turn, half-turn about
    # a cube edge; the starred classes act on displacements as improper operations), typed per class
    fixed = [5, 2, 1, 1, 3]
    trace = [3, 0, -1, -(1 + 0), -(1 - 2)]                      # tr(+-R): 1 + 2 cos theta, negated on the starred classes
    chi_3n = [f * t for f, t in zip(fixed, trace)]
    gamma_3n = {k: exact_mean(sum(n*c*v[i] for i, (n, _, _, _), c in zip(range(5), cls, chi_3n)), 24) for k, v in TD.items()}
    gamma_vib = dict(gamma_3n); gamma_vib['T1'] -= 1; gamma_vib['T2'] -= 1   # minus rotations (T1) and translations (T2)
    assert sum(m * TD[k][0] for k, m in gamma_3n.items()) == 15 and sum(m * TD[k][0] for k, m in gamma_vib.items()) == 9
    # D^J restricted to H through h -> R_h^-1, chi_J = sum_{m=-J}^{J} cos(m theta) by the Chebyshev recursion from the
    # classes' cos theta (turns 0, 120, 180, 90, 180 degrees), exactly; and the physical states per J by parity
    cos_turn = [Fraction(1), Fraction(-1, 2), Fraction(-1), Fraction(0), Fraction(-1)]
    def chi_J(J, c):
        cm = [Fraction(1), c]
        while len(cm) <= J: cm.append(2*c*cm[-1] - cm[-2])
        return cm[0] + 2*sum(cm[1:J + 1])
    j_table = []
    for J in range(0, 7):
        ch = [chi_J(J, c) for c in cos_turn]
        mult = {k: Fraction(sum(n*c*v[i] for i, (n, _, _, _), c in zip(range(5), cls, ch)), 24) for k, v in TD.items()}
        assert all(m.denominator == 1 for m in mult.values()) and sum(m*TD[k][0] for k, m in mult.items()) == 2*J + 1
        mult = {k: int(m) for k, m in mult.items()}
        j_table.append({'J': J, 'species': {k: m for k, m in mult.items() if m},
                        'even': sum(m*weights['Even'][k] for k, m in mult.items()), 'odd': sum(m*weights['Odd'][k] for k, m in mult.items())})
    # under the proper rotations T = A4: classes E, 3 half-turns, 4 and 4 third-turns, with chi_spin 16, 4, 4, 4; the
    # species A, 1E, 2E, T with values in Q(omega)
    w = Q.QOmega(0, 1)
    TT = {'A': [1, 1, 1, 1], '1E': [1, 1, w, w*w], '2E': [1, 1, w*w, w], 'T': [3, -1, 0, 0]}
    tcls = [1, 3, 4, 4]; tspin = [16, 4, 4, 4]
    spin_t = {k: int((sum((Q.QOmega(n*sp)*Q.QOmega(v).conj() for n, sp, v in zip(tcls, tspin, ch)), Q.QOmega(0))*Q.QOmega(Fraction(1, 12))).rational()) for k, ch in TT.items()}
    isomers = [('A', 'A', 1), ('1E', '2E', 1), ('2E', '1E', 1), ('T', 'T', 3)]      # (Gamma_rot, Gamma_nuc, d): Gamma_rot x Gamma_nuc contains A
    kernel_index = {k: 12 // sum(n for n, v in zip(tcls, ch) if Q.QOmega(v) == Q.QOmega(ch[0])) for k, ch in TT.items()}
    summary['methane_rigid'] = {'spin_Td': spin_td, 'weights': weights, 'spin_T': spin_t,
                                'isomers': [(r, nu, d, spin_t[nu]) for r, nu, d in isomers],
                                'entangled_fraction': [3*spin_t['T'], 16], 'monodromy_orders': kernel_index,
                                'gamma_3N': gamma_3n, 'gamma_vib': gamma_vib, 'J_table': j_table}

    # ---- methylamine
    X0 = G.methylamine(0.0, IOTA0)
    emit_pic(mol, 'mol3d-mla-ref', pic_code_rods(X0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, labels=True))
    # figure 5 turns its glyphs in the page so that the projected half-turn axis of b (the body y axis) stands upright
    ax = G.project([(0.0, 1.0, 0.0)], CAM_MLA)[0]
    num.append(f'\\def\\mlaTurn{{{-math.degrees(math.atan2(-ax[0], ax[1])):.2f}}}')
    emit_pic(mol, 'mol3d-water-X', pic_code_rods(wX, G.WATER_ELEMENTS, bonds_w, CAM_FLAT, labels=True))
    emit_pic(mol, 'mol3d-water-RX', pic_code_rods(G.rotate(Rw, wX), G.WATER_ELEMENTS, bonds_w, CAM_FLAT, labels=True))
    emit_pic(mol, 'mol3d-water-minusX', pic_code_rods(G.invert_config(wX), G.WATER_ELEMENTS, bonds_w, CAM_FLAT, labels=True))
    for k, (d, th) in enumerate(samples):
        emit_pic(mol, f'mol3d-water-sample-{k+1}', pic_code_rods(G.water(d, th), G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    # b X0 and R_b X0
    bX0 = G.apply_perm_inversion(X0, Q.b[0], Q.b[1])
    emit_pic(mol, 'mol3d-mla-bref', pic_code_rods(bX0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, labels=True))
    # Newman projections along the t path: tau = 0 (H), pi/3 (eclipsed), 2 pi/3 (tH), and tuH (figure 6)
    for deg in (0, 60, 120):
        emit_pic(mol, f'newman-tau-{deg}', pic_code_newman(G.methylamine(math.radians(deg), IOTA0), twist_deg=(20.0 if deg == 60 else 0.0)))
    emit_pic(mol, 'newman-tuH', pic_code_newman(G.methylamine(-math.pi/3, -IOTA0)))
    # fragment loops in page coordinates (figure 5)
    P0 = G.project(X0, CAM_MLA)
    for name, members in (('methyl', [1, 2, 3, 6]), ('amino', [4, 5, 7])):
        pts = rounded_loop([(P0[i-1][0], P0[i-1][1]) for i in members], 0.78)
        num.append(f'\\def\\fragloop{name}{{' + ' '.join(f'({fmt(x)},{fmt(y)})' for x, y in pts) + '}')
    # version labels (figure 5): position j carries label g(j)
    versions = version_positions()
    vlines = [f'{{{cname}}}/{{{gname}}}/' + '/'.join(str(g[0][j] + 1) for j in range(7)) for cname, gname, g, ta, io in versions]
    num.append('\\def\\versionList{' + ','.join(vlines) + '}')
    summary['versions'] = [(c, gname, round(ta/math.pi, 6), round(io, 5)) for c, gname, g, ta, io in versions]
    # the H-invariant cell U = [-pi/3, pi/3] x [0, iota0] (tau in units of pi, iota in units of iota0) and its translates,
    # the cell of the version gH shifted by its tau and on its side of iota = 0; a cell across the seam is split
    cells = []
    for cname, gname, g, ta, io in versions:
        lo = ((ta/math.pi - 1/3 + 1) % 2) - 1; hi = lo + 2/3
        side = (0, 1) if io > 0 else (-1, 0)
        cells += [(lo, min(hi, 1.0)) + side] + ([(-1.0, hi - 2) + side] if hi > 1 + 1e-9 else [])
    num.append('\\def\\chartCells{' + ','.join(f'{a:.4f}/{b:.4f}/{e1}/{e2}' for a, b, e1, e2 in cells) + '}')

    # normal frame and T_tu (figure 10)
    nf = normal_frame_data()
    summary['normal_frame'] = {'ortho': nf['ortho'], 'M': nf['M'], 'resid': nf['resid'], 'tangent_dim': nf['tangent_dim']}
    amp = 0.42        # the normal displacement drawn (angstrom), illustrative
    CAM_SIDE = G.camera(azimuth_deg=75.0, elevation_deg=18.0)   # the C-N axis lies in the page: both normal displacements at full length;
    # azimuth 75 (not 90) so that no methyl hydrogen sits on the line of sight through the carbon
    def arrows_side(e, longest=1.1):
        # side views: each direction's arrows scaled so that its largest arrow is `longest` angstrom on the page
        # (the direction is the content; a common gain cannot serve a hydrogen twist and a heavy-atom stretch), a drawing choice
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
    emit_pic(mol, 'mla-rot-disp', pic_code_rods(G.rotate(Rdisp, Xdisp), G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))

    # generic X and bX (figure 11): X = X0(0.22, iota0) + a small normal displacement
    Xg = G.methylamine(0.22, IOTA0)
    eg, _ = G.normal_vector([(0,0,0)]*3 + [(0,0,1.0), (0,0,-1.0)] + [(0,0,0)]*2, 0.22, IOTA0)
    Xg = [G.add(x, G.scale(v, 0.30)) for x, v in zip(Xg, eg)]
    bXg = G.apply_perm_inversion(Xg, Q.b[0], Q.b[1])
    R, res = G.kabsch_rotation(Xg, bXg, G.MLA_MASSES)
    summary['generic_X_bX_best_rotation_residual'] = res
    emit_pic(mol, 'mol3d-mla-X', pic_code_rods(Xg, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, labels=True))
    emit_pic(mol, 'mol3d-mla-bX', pic_code_rods(bXg, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, labels=True))

    # ---- KRb (figure 4)
    K = G.krb_schematic()
    Sk, Gin, GK, GRb, _, _, Es = Q.krb_groups()
    # the three channel groups meet exactly in <E*> and cover S
    core = Gin & GK & GRb
    assert core == (Gin & GK) == (Gin & GRb) == (GK & GRb) == {Q.E(4), Es}, 'channel groups must meet in <E*>'
    assert Gin | GK | GRb == Sk, 'channel groups must cover S'
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
    # the named subgroups at the center of their rows: G_Rb in the middle, G_in at its left, G_K at its right, so that
    # the strands run up from G_in to S and down to G_K and G_Rb without crossing
    priority = {frozenset(GRb): 0, frozenset(Gin): 1, frozenset(GK): 2}
    for L, row in rows.items():
        row.sort(key=lambda K_: (0 if K_ in named else 1, priority.get(K_, 0), sorted(K_)))
        n_ = len(row); named_in = [K_ for K_ in row if K_ in named]; others = [K_ for K_ in row if K_ not in named]
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
            nodes.append((key_[K_], x, ylev[L], 1 if K_ in named else 0))
    cov_ = []
    for A in allsubs:
        for B in allsubs:
            if A < B and len(B) == 2*len(A):
                cov_.append((key_[A], key_[B], 1 if (A in named and B in named) else 0))
    num.append('\\def\\krbNodes{' + ','.join(f'{i}/{x:.4f}/{y:.1f}/{m}' for i, x, y, m in nodes) + '}')
    idx = {named[K_]: key_[K_] for K_ in named}
    core_name = '\\langle E^*\\rangle'; in_name = 'G_{\\mathrm{in}}'; k_name = 'G_{\\mathrm K}'; rb_name = 'G_{\\mathrm{Rb}}'
    num.append('\\def\\krbIdxE{%d}\\def\\krbIdxCore{%d}\\def\\krbIdxIn{%d}\\def\\krbIdxK{%d}\\def\\krbIdxRb{%d}\\def\\krbIdxS{%d}'
               % (idx['E'], idx[core_name], idx[in_name], idx[k_name], idx[rb_name], idx['S']))
    num.append('\\def\\krbCovers{' + ','.join(f'{a}/{b}/{m}' for a, b, m in cov_) + '}')

    # ---- the subgroup intervals (figure 5)
    subs, names, cov = Q.interval_data()
    order_of = {K: len(K)//2 for K in subs}     # |K/H|
    # drawing x-slots per node (design choice); y = log2 of the version count |K/H|
    xslot = {'H': 0.0, '\\langle H,u\\rangle': 4.6, '\\langle H,E^*\\rangle': -0.4, '\\langle H,u^*\\rangle': -2.9,
             'G_6': 2.6, '\\langle H,u,E^*\\rangle': -2.6, 'G_{12}': 3.2, '\\langle G_6,E^*\\rangle': -0.6,
             '\\langle G_6,u^*\\rangle': -4.4, 'B': 0.0}
    ordered = sorted(subs, key=lambda K: (len(K), names[K]))
    key = {names[K]: i for i, K in enumerate(ordered)}
    num.append('\\def\\hasseNodes{' + ','.join(f'{key[names[K]]}/{xslot[names[K]]}/{math.log2(len(K)//2):.4f}' for K in ordered) + '}')
    # the generator added along each cover of [H, B], and a minimal generator form of every subgroup (b always first)
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
                for combo in itertools.combinations(cands, r):
                    if frozenset(Q.close([Q.b] + [g for _, g in combo], Q.N)) == frozenset(K_):
                        found = ','.join(['b'] + [nm for nm, _ in combo]); break
                if found: break
        assert found, 'no generator form'
        gens[key[names[K_]]] = found
    labs = []
    for K_ in subs:
        k = key[names[K_]]; gform = '\\langle ' + gens[k] + '\\rangle'
        labs.append((k, (names[K_] + '=' + gform) if names[K_] in chain_forms else gform))
    num.append('\\def\\hasseLabels{' + ','.join(f'{k}/{{{l}}}' for k, l in sorted(labs)) + '}')
    # the whole interval [H, S]: every subgroup containing H, with covers; the bond interval is a sub-poset
    big = Q.interval(Q.H, Q.S, Q.N)
    bigcov = Q.covers(big)
    bigsorted = sorted(big, key=lambda K: (len(K), sorted(K)))
    bigkey = {K: i for i, K in enumerate(bigsorted)}
    inbond = {bigkey[K] for K in big if K <= Q.B}
    levels = sorted(set(len(K) for K in big))
    yof = {K: math.log(len(K)/2)/math.log(120) for K in big}      # height by log of the version count
    below = {K: [a for a, b in bigcov if b == K] for K in big}
    above = {K: [b for a, b in bigcov if a == K] for K in big}
    rows = {L: [K for K in bigsorted if len(K) == L] for L in levels}

    def layout(pull, median, sweeps):
        """Layered placement: y by level; x by repeated barycenter (or median) ordering of each level from its
        neighbors, alternating the sweep direction, nodes spread evenly in each level; the bond interval pulled
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
    # four nodes nudged so that no cover passes through a node: G_6 and another node of the bond interval, two outside it
    for k_, x_ in ((6, 0.38), (8, 0.12), (15, 0.21), (26, 0.26)):
        xpos[bigsorted[k_]] = x_
    assert names.get(bigsorted[6]) == 'G_6' and bigsorted[8] <= Q.B and not bigsorted[15] <= Q.B and not bigsorted[26] <= Q.B
    num.append('\\def\\bigNodes{' + ','.join(f'{bigkey[K]}/{xpos[K]:.4f}/{yof[K]:.4f}/{1 if bigkey[K] in inbond else 0}' for K in bigsorted) + '}')
    num.append('\\def\\bigCovers{' + ','.join(f'{bigkey[a]}/{bigkey[b]}/{1 if (a <= Q.B and b <= Q.B) else 0}' for a, b in bigcov) + '}')
    chain_keys = {names[K]: bigkey[K] for K in big if K in names and names[K] in ('H', 'G_6', 'G_{12}', 'B')}
    num.append('\\def\\bigChain{' + ','.join(str(chain_keys[n]) for n in ('H', 'G_6', 'G_{12}', 'B')) + '}')
    chain_seq = [chain_keys[n] for n in ('H', 'G_6', 'G_{12}', 'B')]
    num.append('\\def\\bigChainCovers{' + ','.join(f'{a}/{b}' for a, b in zip(chain_seq, chain_seq[1:])) + '}')
    summary['interval_HS'] = {'subgroups': len(big), 'covers': len(bigcov), 'in_bond_interval': len(inbond)}
    summary['interval'] = {'subgroups': len(subs), 'covers': len(cov),
                           'versions': sorted(order_of.values())}

    # ---- spin weights (figure 8)
    w = Q.spin_weights()
    summary['spin_weights'] = {str(k): v for k, v in w.items()}
    num.append(f'\\def\\spinAone{{{w[1]["A1"]}}}\\def\\spinAtwo{{{w[1]["A2"]}}}\\def\\spinE{{{w[1]["E"]}}}\\def\\spinDim{{{w[1]["A1"]+w[1]["A2"]+2*w[1]["E"]}}}')   # figure 8, the proton factor

    # ---- Gaussian shares (figure 9)
    with open(os.path.join(DATA, 'gaussian-shares.dat'), 'w') as f:
        f.write('x wA1 wE\n')
        f.write('0.00 0.333333 0.666667\n')   # the limit: no overlap, shares 1/3 and 2/3
        for i in range(2, 101):
            x = i/100
            c = math.exp(-1/(8*x*x))
            f.write(f'{x:.2f} {(1+2*c)/3:.6f} {2*(1-c)/3:.6f}\n')
    # figure 9's three rows: Delta/d, the packet width in degrees for d = 120 degrees, the A1 share
    rows = [(r, (1 + 2*math.exp(-1/(8*float(r)**2)))/3) for r in (Fraction(1, 6), Fraction(1, 3), Fraction(2, 3))]
    num.append('\\def\\packetRows{' + ','.join(f'{{{r}}}/{float(r)*120:g}/{w:.4f}' for r, w in rows) + '}')
    num.append('\\def\\packetMarks{' + ' '.join(f'({float(r):.4f},{w:.4f})' for r, w in rows) + '}')

    # ---- component functions (figure 2): a physical choice on the delta = 0 slice.  There X P_sigma = R X with R the
    # in-plane half-turn, so a rotation-invariant state obeys Psi(X) = chi_stat(sigma) sigma.Psi(X), i.e. f_zeta = -f_xi:
    # the proton singlet (xi - zeta) times a scalar f(theta), f an illustrative Gaussian in theta.
    COMP_CENTER, COMP_WIDTH = 104.5, 14.0   # degrees
    num.append(f'\\def\\compCenter{{{COMP_CENTER}}}\\def\\compWidth{{{COMP_WIDTH}}}')
    def f_xi(th): return math.exp(-((th-COMP_CENTER)/COMP_WIDTH)**2/2)
    def f_zeta(th): return -f_xi(th)   # the singlet: f_zeta = -f_xi
    with open(os.path.join(DATA, 'component-functions.dat'), 'w') as f:
        f.write('theta fxi fzeta norm\n')
        for i in range(0, 101):
            th = 80 + i
            f.write(f'{th} {f_xi(th):.5f} {f_zeta(th):.5f} {f_xi(th)**2+f_zeta(th)**2:.5f}\n')
    for k, (d, th) in enumerate(samples):
        tag = 'one two three'.split()[k]
        num.append(f'\\def\\sampleTheta{tag}{{{th:g}}}')

    # ---- momentum labels (figure 12)
    rho = Fraction(1, 2)   # the symmetric frame: each turn of the frame is half the torsion (Mellor, Yurchenko, Mant, Jensen 2019)
    num.append(f'\\def\\rhoIll{{{float(rho)}}}')
    rows = []
    for Kq in (-2, -1, 0, 1, 2):
        rows.append(f'{Kq}/' + '{' + ','.join(f'{float(m + rho*Kq):.1f}' for m in range(-3, 4)) + '}')
    num.append('\\def\\kappaRows{' + ','.join(rows) + '}')
    tms = torsion_mode_species()
    num.append('\\def\\torsionSpecies{' + ','.join(f'{m}/{c}/{s if s else "none"}' for m, c, s in tms) + '}')
    summary['torsion_species'] = tms

    with open(os.path.join(DATA, 'molecules.tex'), 'w') as f:
        f.write('\n'.join(mol) + '\n')
    with open(os.path.join(DATA, 'numbers.tex'), 'w') as f:
        f.write('\n'.join(num) + '\n')
    with open(os.path.join(DATA, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=1, default=float)
    print('wrote', DATA)

if __name__ == '__main__':
    main()
