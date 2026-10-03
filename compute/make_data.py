"""Generate the verified figure data in figures/data/ from compute/.

Run from the repository root:  python3 compute/make_data.py
Everything emitted here is recomputed from declared masses, model parameters
and group definitions; nothing is typed in by hand.  checks/verify.py
re-derives the same quantities independently and asserts them.
"""
import json, math, os, sys
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

def pic_code(X, elements, bonds, cam, arrows=None, label_macros=None):
    """TikZ pic code: depth-sorted half-bonds and atom discs, in angstrom units."""
    P = G.project(X, cam)
    order = G.draw_order(P)
    n = len(X)
    if label_macros is None:
        label_macros = ['\\mol' + 'ABCDEFG'[i] for i in range(n)]
    nbrs = {i: [] for i in range(n)}
    for a, b in bonds:
        nbrs[a-1].append(b-1); nbrs[b-1].append(a-1)
    lines = []
    for i in order:
        px, py, _ = P[i]
        for j in nbrs[i]:
            qx, qy, _ = P[j]
            mx, my = 0.5*(px+qx), 0.5*(py+qy)
            lines.append(f'\\molbond{{{fmt(px)}}}{{{fmt(py)}}}{{{fmt(mx)}}}{{{fmt(my)}}}')
        lines.append(f'\\molatom{{{elements[i]}}}{{{label_macros[i]}}}{{{fmt(px)}}}{{{fmt(py)}}}')
    # displacement arrows go on top of every disc so that none is hidden
    if arrows is not None:
        for i in order:
            if arrows[i] is not None:
                px, py, _ = P[i]; dx, dy = arrows[i]
                lines.append(f'\\molarrow{{{fmt(px)}}}{{{fmt(py)}}}{{{fmt(dx)}}}{{{fmt(dy)}}}')
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

    The candidate formulas are the ones derived in docs/verification.md; here
    they are confirmed numerically at several (tau, eta) by a residual check.
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
    emit_pic(mol, 'water-X', pic_code(wX, G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    Rw = G.rot_z(math.radians(40.0))
    emit_pic(mol, 'water-RX', pic_code(G.rotate(Rw, wX), G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    emit_pic(mol, 'water-minusX', pic_code(G.invert_config(wX), G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    for i, c in enumerate(wX):
        num.append(f'\\def\\waterX{"ABC"[i]}x{{{fmt(c[0],2)}}}\\def\\waterX{"ABC"[i]}y{{{fmt(c[1],2)}}}')
    num.append(f'\\def\\massH{{{G.MASS["H"]:.3f}}}\\def\\massO{{{G.MASS["O"]:.3f}}}')
    num.append(f'\\def\\waterRotDeg{{40}}')
    # shape grid for figure 1(c): delta ratios and angles, drawn mass-centred
    grid_d = [-0.5, -0.25, 0.0, 0.25, 0.5]; grid_t = [60.0, 90.0, 120.0, 150.0, 180.0]
    for i, d in enumerate(grid_d):
        for j, th in enumerate(grid_t):
            Xg = G.water(d, th)
            emit_pic(mol, f'water-grid-{i}{j}', pic_code(Xg, G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    summary['water_grid'] = {'delta_ratios': grid_d, 'theta_deg': grid_t}
    # three sample configurations for figure 2
    samples = [(0.0, 104.5), (0.15, 125.0), (-0.1, 95.0)]
    for k, (d, th) in enumerate(samples):
        emit_pic(mol, f'water-sample-{k+1}', pic_code(G.water(d, th), G.WATER_ELEMENTS, bonds_w, CAM_FLAT))
    summary['water_samples'] = samples

    # ---- methylamine
    X0 = G.methylamine(0.0, ETA0)
    summary['methylamine_frame'] = {k: round(v, 5) for k, v in G.MLA_FRAME.items()}
    summary['methylamine_X0'] = [[round(v, 4) for v in c] for c in X0]
    summary['methylamine_X0_Xm'] = [round(v, 12) for v in G.mass_moment(X0, G.MLA_MASSES)]
    emit_pic(mol, 'mla-ref', pic_code(X0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    # b X0 and R_b X0
    bX0 = G.apply_perm_inversion(X0, Q.b[0], Q.b[1])
    Rb = G.rot_y(math.pi)
    res_b = G.config_distance(bX0, G.rotate(Rb, X0))
    summary['b_equals_rotation_residual'] = res_b
    emit_pic(mol, 'mla-bref', pic_code(bX0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    # rigid orbit samples: the same shape X0 in two other orientations (figure 5a)
    for tag, axis, ang in (('a', (0.2, 1.0, 0.3), 70.0), ('b', (1.0, 0.1, -0.4), 150.0)):
        emit_pic(mol, f'mla-rot-{tag}', pic_code(G.rotate(G.rot_axis(axis, math.radians(ang)), X0), G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    # shapes along the t path at fixed orientation: tau = pi/3 (eclipsed) and 2 pi/3 (the version tH) (figure 6d)
    for deg in (60, 120):
        emit_pic(mol, f'mla-tau-{deg}', pic_code(G.methylamine(math.radians(deg), ETA0), G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    summary['tau_shapes_deg'] = [0, 60, 120]
    # fragment loops in page coordinates (figure 4)
    P0 = G.project(X0, CAM_MLA)
    for name, members in (('methyl', [1, 2, 3, 6]), ('amino', [4, 5, 7]), ('all', [1, 2, 3, 4, 5, 6, 7])):
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
        num.append(f'\\def\\versionTau{gname.replace("^2","sq")}{{{fmt(ta/math.pi,4)}}}')
    num.append('\\def\\versionList{' + ','.join(vlines) + '}')
    summary['versions'] = [(c, g, round(ta/math.pi, 6), round(et, 5)) for c, g, ta, et in versions]
    actions, worst = derive_family_actions()
    summary['family_actions'] = actions; summary['family_action_residual'] = worst
    num.append(f'\\def\\etaZero{{{fmt(ETA0,2)}}}')

    # normal frame and T_tu (figure 10)
    nf = normal_frame_data()
    summary['normal_frame'] = {'gram': nf['gram'], 'ortho': nf['ortho'], 'M': nf['M'],
                               'resid': nf['resid'], 'cfg_res': nf['cfg_res'],
                               'tangent_dim': nf['tangent_dim']}
    amp = 0.42        # actual displacement amplitude Q (angstrom) used for the displaced drawing
    arrow_gain = 3.0  # arrows are drawn 2.5 x longer than the displacement; declared in the figure
    def arrows_for(e):
        Pe = [(G.dot(v, CAM_MLA[0]), G.dot(v, CAM_MLA[1])) for v in e]
        return [(arrow_gain*amp*px, arrow_gain*amp*py) if math.hypot(px, py)*amp > 0.12 else None for px, py in Pe]
    emit_pic(mol, 'mla-ref-eplus', pic_code(X0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, arrows=arrows_for(nf['e_plus'])))
    emit_pic(mol, 'mla-ref-eminus', pic_code(X0, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA, arrows=arrows_for(nf['e_minus'])))
    Xdisp = [G.add(x, G.scale(v, amp)) for x, v in zip(X0, nf['e_minus'])]
    Rdisp = G.rot_axis((0.3, 1.0, 0.25), math.radians(55.0))
    emit_pic(mol, 'mla-disp', pic_code(Xdisp, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    emit_pic(mol, 'mla-rot-disp', pic_code(G.rotate(Rdisp, Xdisp), G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    num.append(f'\\def\\dispAmp{{{fmt(amp,2)}}}\\def\\arrowGain{{{arrow_gain:g}}}')

    # generic X and bX (figure 11b): X = X0(0.22, eta0) + small normal displacement
    Xg = G.methylamine(0.22, ETA0)
    eg, _ = G.normal_vector([(0,0,0)]*3 + [(0,0,1.0), (0,0,-1.0)] + [(0,0,0)]*2, 0.22, ETA0)
    Xg = [G.add(x, G.scale(v, 0.30)) for x, v in zip(Xg, eg)]
    bXg = G.apply_perm_inversion(Xg, Q.b[0], Q.b[1])
    R, res = G.kabsch_rotation(Xg, bXg, G.MLA_MASSES)
    summary['generic_X_bX_best_rotation_residual'] = res
    summary['generic_X_Xm'] = [round(v, 12) for v in G.mass_moment(Xg, G.MLA_MASSES)]
    emit_pic(mol, 'mla-X', pic_code(Xg, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    emit_pic(mol, 'mla-bX', pic_code(bXg, G.MLA_ELEMENTS, G.MLA_BONDS, CAM_MLA))
    # free-action check near the reference: min_g ||gX0 - X0|| over g != E in G12
    dmin = min(G.config_distance(G.apply_perm_inversion(X0, g[0], g[1]), X0) for g in Q.G12 if g != Q.E(Q.N))
    summary['min_displacement_by_nontrivial_G12_element'] = dmin
    num.append(f'\\def\\minDispG{{{dmin:.2f}}}\\def\\bXresidual{{{res:.2f}}}')

    # sheets (figure 11a): cosets of G12/H with representatives
    sheet_pairs = [('E', 'b'), ('t', 'tb'), ('t^2', 't^2b'), ('u', 'ub'), ('tu', 'tub'), ('t^2u', 't^2ub')]
    num.append('\\def\\sheetPairs{' + ','.join(f'{{{a}}}/{{{b}}}' for a, b in sheet_pairs) + '}')
    # continuation results
    tb = Q.mul(Q.t, Q.b); bt = Q.mul(Q.b, Q.t); t2b = Q.mul(Q.mul(Q.t, Q.t), Q.b)
    summary['bt_equals_t2b'] = (bt == t2b); summary['tb_ne_bt'] = (tb != bt)
    summary['btb_equals_tinv'] = (Q.mul(Q.b, Q.mul(Q.t, Q.b)) == Q.inv(Q.t))

    # ---- KRb (figure 4b)
    K = G.krb_schematic()
    emit_pic(mol, 'krb-config', pic_code(K, G.KRB_ELEMENTS, [], CAM_FLAT))
    PK = G.project(K, CAM_FLAT)
    for name, parts in (('in', [[1, 3], [2, 4]]), ('complex', [[1, 2, 3, 4]]), ('out', [[1, 2], [3, 4]])):
        for k, members in enumerate(parts):
            pts = rounded_loop([(PK[i-1][0], PK[i-1][1]) for i in members], 1.05)
            num.append(f'\\def\\krbloop{name}{"AB"[k]}{{' + ' '.join(f'({fmt(x)},{fmt(y)})' for x, y in pts) + '}')
    Sk, Gin, GK, GRb, *_ = Q.krb_groups()
    summary['krb_orders'] = [len(Sk), len(Gin), len(GK), len(GRb)]
    summary['krb_any_two_generate'] = all(len(Q.close(list(a | b), 4)) == 8 for a, b in [(Gin, GK), (Gin, GRb), (GK, GRb)])

    # ---- subgroup interval (figure 5b)
    subs, names, cov = Q.interval_data()
    order_of = {K: len(K)//2 for K in subs}     # |K/H|
    # drawing x-slots per node (design choice); y = log2 of the version count |K/H|
    xslot = {'H': 0.0, '\\langle H,u\\rangle': 4.6, '\\langle H,E^*\\rangle': -0.4, '\\langle H,u^*\\rangle': -2.9,
             'G_6': 2.6, '\\langle H,u,E^*\\rangle': -2.6, 'G_{12}': 2.4, '\\langle G_6,E^*\\rangle': 0.5,
             '\\langle G_6,u^*\\rangle': -4.6, 'B': 0.0}
    ordered = sorted(subs, key=lambda K: (len(K), names[K]))
    key = {names[K]: i for i, K in enumerate(ordered)}
    num.append('\\def\\hasseNodes{' + ','.join(
        f'{key[names[K]]}/{xslot[names[K]]}/{math.log2(len(K)//2):.4f}/{len(K)//2}/{{{names[K]}}}' for K in ordered) + '}')
    num.append('\\def\\hasseCovers{' + ','.join(f'{key[names[a]]}/{key[names[b]]}' for a, b in cov) + '}')
    num.append('\\def\\hasseChain{' + ','.join(str(key[n]) for n in ('H', 'G_6', 'G_{12}', 'B')) + '}')
    summary['interval'] = {'subgroups': len(subs), 'covers': len(cov),
                           'versions': sorted(order_of.values())}

    # ---- species and spin weights (figures 7, 8)
    cnames, act_t, act_b = Q.coset_action_table()
    summary['coset_action'] = {'t': act_t, 'b': act_b}
    w = Q.spin_weights()
    summary['spin_weights'] = {str(k): v for k, v in w.items()}
    # multiplicities of species in C[G6/H]: permutation character (3, 0, 1)
    perm_char = {'E': 3, 't': 0, 'b': 1}
    local = {n: int(round(sum(Q.G6_CHARS[n][c]*perm_char[c]*{'E': 1, 't': 2, 'b': 3}[c] for c in 'Etb')/6)) for n in Q.G6_CHARS}
    summary['local_copies'] = local
    # TeX control sequences cannot contain digits: A1 -> Aone, A2 -> Atwo
    spname = {'A1': 'Aone', 'A2': 'Atwo', 'E': 'E'}
    for par, tag in ((1, 'Even'), (-1, 'Odd')):
        for sp in ('A1', 'A2', 'E'):
            num.append(f'\\def\\weight{tag}{spname[sp]}{{{w[par][sp]}}}')
    for sp in ('A1', 'A2', 'E'):
        num.append(f'\\def\\localCopies{spname[sp]}{{{local[sp]}}}')
    num.append(f'\\def\\spinDim{{{2**5}}}')

    # ---- Gaussian shares (figure 9)
    with open(os.path.join(DATA, 'gaussian-shares.dat'), 'w') as f:
        f.write('x wA1 wE wA2\n')
        for i in range(2, 201):
            x = i/100
            c = math.exp(-1/(8*x*x))
            f.write(f'{x:.2f} {(1+2*c)/3:.6f} {2*(1-c)/3:.6f} 0\n')
    x_show = 1/3; c_show = math.exp(-1/(8*x_show*x_show))
    num.append(f'\\def\\shareShownA{{{(1+2*c_show)/3:.3f}}}\\def\\shareShownE{{{2*(1-c_show)/3:.3f}}}')
    summary['gaussian_shown'] = {'sigma_over_d': x_show, 'c': c_show, 'wA1': (1+2*c_show)/3, 'wE': 2*(1-c_show)/3}

    # ---- illustrative component functions (figure 2b): values of a continuous representative
    def comp_a(th): return math.exp(-((th-104.5)/14.0)**2/2)
    def comp_b(th): return 0.9*((th-118.0)/16.0)*math.exp(-((th-118.0)/16.0)**2/2)
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
    rho = 0.3
    num.append(f'\\def\\rhoIll{{{rho}}}')
    rows = []
    for Kq in (-1, 0, 1):
        rows.append(f'{Kq}/' + '{' + ','.join(f'{m + rho*Kq:.1f}' for m in range(-3, 4)) + '}')
    num.append('\\def\\kappaRows{' + ','.join(rows) + '}')
    tms = torsion_mode_species()
    num.append('\\def\\torsionSpecies{' + ','.join(f'{m}/{c}/{s if s else "none"}' for m, c, s in tms) + '}')
    summary['torsion_species'] = tms
    summary['kappa_rule'] = 'kappa = m + rho K from (r,tau+2pi,Q) ~ (r R_z(-2 pi rho),tau,Q) and D^J_MK(r R_z(a)) = e^{-iKa} D^J_MK(r)'

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
