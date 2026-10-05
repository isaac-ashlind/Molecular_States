"""The closing pages: rigid methane run through the pipeline and set against Albert, Kubischta, Lemeshko and Liu
(arXiv:2403.04572v4), recomputed from explicit matrices and compared with the numbers make_data emits
(summary.json, keys methane_rigid and methane_loop).

  geometry   every element of Td(M) = {even relabellings} u {odd relabellings}* carries the tetrahedron X0 into its
             rotational orbit (a proper rotation R_h, residual ~ 0); all 24 close as the cube group O, the 12 even as T
  spin       (C^2)^4 under the relabellings: 5 A1 + E + 3 T2 (Td labels); under the proper rotations: 5 A + 1E + 2E + 3 T
  pairing    physical states per rotational species and parity (chi_stat chi_pm is A2 for even parity, A1 for odd)
  vibrations the fifteen Cartesian displacements: A1 + E + T1 + 3 T2; the nine normal to the orbit in the centred
             space: A1 + E + 2 T2 (from explicit 15 x 15 matrices and a mass-orthonormal normal basis)
  J ladder   D^J restricted to H through h -> R_h^-1 (from the actual turning angles), J = 0..6, and the physical
             states of each J by parity
  ball       the Voronoi cell of the identity among the twelve R_h (bi-invariant metric) is the octahedron
             |x|+|y|+|z| <= 1 in Rodrigues coordinates; its vertices are the quarter-turns (the starred 4-cycles);
             the face toward R_g is glued to the opposite face by R -> R R_g^-1 with a third of a twist
  packets    the projection shares of a packet with class overlaps c3, c2 against the closed forms
  isomers    Gamma_rot x Gamma_nuc contains A: (A, A) 5, (1E, 2E) 1, (2E, 1E) 1, (T, T) 3, with d = 1, 1, 1, 3
  entangled  the invariant vector of T x T (T = the rotation matrices R_h) has Schmidt rank 3; 9 of the 16 spin states
  monodromy  |T / ker Gamma| = 1, 3, 3, 12; the loop about the bond to hydrogen 1 returns X0 with 2, 3, 4 cycled
Requires numpy. Exit status nonzero on any mismatch.
"""
import itertools, json, os, sys, cmath, math
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'compute'))
import geometry as G

def parity(p):
    return (-1) ** sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])

def cycle_type(p):
    seen, lens = set(), []
    for i in range(4):
        if i in seen: continue
        n, j = 0, i
        while j not in seen:
            seen.add(j); j = p[j]; n += 1
        lens.append(n)
    return tuple(sorted(lens, reverse=True))

def spin_matrix(p):
    M = np.zeros((16, 16))
    for idx in range(16):
        bits = [(idx >> k) & 1 for k in range(4)]
        new = [0] * 4
        for i in range(4): new[p[i]] = bits[i]
        M[sum(b << k for k, b in enumerate(new)), idx] = 1
    return M

def axis_angle(Rm):
    ang = math.acos(max(-1.0, min(1.0, (np.trace(Rm) - 1) / 2)))
    if ang < 1e-9: return np.zeros(3), 0.0
    if abs(ang - math.pi) < 1e-9:
        M = (Rm + np.eye(3)) / 2; k = int(np.argmax(np.diag(M))); return M[:, k] / math.sqrt(M[k, k]), ang
    n = np.array([Rm[2, 1] - Rm[1, 2], Rm[0, 2] - Rm[2, 0], Rm[1, 0] - Rm[0, 1]]) / (2 * math.sin(ang))
    return n, ang

def main():
    ok = True
    X0 = G.methane_c2()
    masses = np.array(G.CH4_MASSES)
    perms = [list(p) for p in itertools.permutations(range(4))]
    # the stabilizer: sigma unstarred when even, starred when odd; R_h from an exact orthogonal Procrustes fit by SVD
    # (the suite's Newton polar routine is tuned to small turns and half-turns and does not converge for third-turns)
    def fit(X, Y):
        A = np.array(X).T; B = np.array(Y).T; W = np.diag(G.CH4_MASSES)
        U, _, Vt = np.linalg.svd(B @ W @ A.T)
        D = np.diag([1, 1, np.sign(np.linalg.det(U @ Vt))])
        Rm = U @ D @ Vt
        return Rm, float(np.sqrt(np.sum(W @ ((Rm @ A - B) ** 2).T)))
    R = {}
    for p in perms:
        star = parity(p) < 0
        Y = G.apply_perm_inversion(X0, p + [4], star)
        Rh, res = fit(X0, Y)
        ok &= res < 1e-9 and abs(np.linalg.det(Rh) - 1) < 1e-9
        R[tuple(p)] = Rh
    even = [tuple(p) for p in perms if parity(p) > 0]
    allp = [tuple(p) for p in perms]
    closed = all(any(np.allclose(R[a] @ R[b], R[c]) for c in even) for a in even for b in even)
    closed_all = all(any(np.allclose(R[a] @ R[b], R[c]) for c in allp) for a in allp for b in allp)
    angles = sorted(round(math.degrees(axis_angle(R[q])[1])) for q in allp)
    ok &= closed and len(even) == 12 and closed_all and angles == [0] + [90] * 6 + [120] * 8 + [180] * 9
    print('geometry: 24 relabellings carry X0 into its orbit by proper rotations; the 12 even ones close:', closed,
          '| all 24 close as the cube group O (1, 6 quarter-turns, 8 third-turns, 9 half-turns):', closed_all and angles == [0] + [90] * 6 + [120] * 8 + [180] * 9)
    # spin characters and the Td decomposition (A1, A2, E, T1, T2 by cycle type and star)
    TD = {'A1': {(1,1,1,1): 1, (3,1): 1, (2,2): 1, (4,): 1, (2,1,1): 1},
          'A2': {(1,1,1,1): 1, (3,1): 1, (2,2): 1, (4,): -1, (2,1,1): -1},
          'E':  {(1,1,1,1): 2, (3,1): -1, (2,2): 2, (4,): 0, (2,1,1): 0},
          'T1': {(1,1,1,1): 3, (3,1): 0, (2,2): -1, (4,): 1, (2,1,1): -1},
          'T2': {(1,1,1,1): 3, (3,1): 0, (2,2): -1, (4,): -1, (2,1,1): 1}}
    chi = {tuple(p): np.trace(spin_matrix(p)) for p in perms}
    spin_td = {k: round(sum(chi[q] * v[cycle_type(q)] for q in chi) / 24) for k, v in TD.items()}
    ok &= spin_td == {'A1': 5, 'A2': 0, 'E': 1, 'T1': 0, 'T2': 3}
    print('spin under Td(M):', spin_td)
    # pairing by parity: chi_total(h) = sign(h) * (parity if starred)
    weights = {}
    for par, tag in ((1, 'Even'), (-1, 'Odd')):
        weights[tag] = {k: round(sum(v[cycle_type(q)] * chi[q] * parity(q) * (par if parity(q) < 0 else 1) for q in chi) / 24) for k, v in TD.items()}
    ok &= weights['Even'] == {'A1': 0, 'A2': 5, 'E': 1, 'T1': 3, 'T2': 0} and weights['Odd'] == {'A1': 5, 'A2': 0, 'E': 1, 'T1': 0, 'T2': 3}
    print('physical states per rotational species, even | odd parity:', weights)
    # vibrations: h acts on a displacement dX (3 x 5) by dX -> +-R_h^-1 dX P_sigma (relabel, invert if starred, rotate
    # back into the body frame); explicit 15 x 15 matrices; then restricted to the centred space minus the orbit tangent
    def disp_matrix(q):
        P = np.array(G.perm_matrix(list(q) + [4], 5), float)        # X P : column j <- column of the relabelling
        sgn = -1.0 if parity(q) < 0 else 1.0
        Rt = R[q].T
        M = np.zeros((15, 15))
        for i in range(5):
            for a in range(3):
                E = np.zeros((3, 5)); E[a, i] = 1.0
                M[:, 3 * i + a] = (sgn * Rt @ E @ P).T.reshape(15)
        return M
    chi_3n = {q: np.trace(disp_matrix(q)) for q in allp}
    gamma_3n = {k: round(sum(chi_3n[q] * v[cycle_type(q)] for q in allp) / 24) for k, v in TD.items()}
    ok &= gamma_3n == {'A1': 1, 'A2': 0, 'E': 1, 'T1': 1, 'T2': 3}
    # the normal space: mass-orthogonal to the three translations and the three rotation generators
    Wm = np.diag(np.repeat(masses, 3))
    trans = [np.tile(np.eye(3)[a], 5) for a in range(3)]
    rots = [np.array(g).reshape(15) for g in G.rotation_generators(X0)]
    span = np.array(trans + rots).T
    # mass-orthonormal basis of the complement: solve in the metric Wm
    Q, _ = np.linalg.qr(np.sqrt(Wm) @ span)
    comp = np.eye(15) - Q @ Q.T                                   # projector in the mass-weighted coordinates
    U, s, _ = np.linalg.svd(comp)
    N = U[:, :9] if np.sum(s > 1e-9) == 9 else None
    ok &= N is not None
    inv_sqrt = np.diag(1 / np.sqrt(np.repeat(masses, 3)))
    chi_vib = {q: np.trace(N.T @ np.sqrt(Wm) @ disp_matrix(q) @ inv_sqrt @ N) for q in allp}
    gamma_vib = {k: round(sum(chi_vib[q] * v[cycle_type(q)] for q in allp) / 24) for k, v in TD.items()}
    ok &= gamma_vib == {'A1': 1, 'A2': 0, 'E': 1, 'T1': 0, 'T2': 2}
    print('displacements under Td(M):', gamma_3n, '| the nine normal to the orbit:', gamma_vib)
    # the J ladder: D^J restricted to H through h -> R_h^-1 (same angle as R_h), by characters
    def chi_J(J, Rm):
        ang = axis_angle(Rm)[1]
        return 2 * J + 1 if ang < 1e-9 else math.sin((2 * J + 1) * ang / 2) / math.sin(ang / 2)
    j_table = []
    for J in range(7):
        mult = {k: round(sum(chi_J(J, R[q]) * v[cycle_type(q)] for q in allp) / 24) for k, v in TD.items()}
        ok &= sum(m * TD[k][(1,1,1,1)] for k, m in mult.items()) == 2 * J + 1
        j_table.append({'J': J, 'species': {k: m for k, m in mult.items() if m},
                        'even': sum(m * weights['Even'][k] for k, m in mult.items()), 'odd': sum(m * weights['Odd'][k] for k, m in mult.items())})
    ok &= j_table[0]['species'] == {'A1': 1} and j_table[1]['species'] == {'T1': 1} and j_table[2]['species'] == {'E': 1, 'T2': 1}
    ok &= j_table[3]['species'] == {'A2': 1, 'T1': 1, 'T2': 1} and j_table[4]['species'] == {'A1': 1, 'E': 1, 'T1': 1, 'T2': 1}
    print('J ladder:', [(r['J'], r['species'], r['even'], r['odd']) for r in j_table[:5]])
    # the orientation ball: Voronoi cell of the identity among the twelve R_h = the octahedron in Rodrigues coordinates
    rng = np.random.default_rng(11)
    def rodrigues(Rm):
        n, ang = axis_angle(Rm); return math.tan(ang / 2) * n
    def angle_of(Rm): return axis_angle(Rm)[1]
    bad = 0; tested = 0
    for _ in range(3000):
        qv = rng.normal(size=4); qv /= np.linalg.norm(qv); w, x, y, z = qv
        Rm = np.array([[1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)], [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)], [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)]])
        rod = rodrigues(Rm); l1 = np.abs(rod).sum()
        if abs(l1 - 1) < 2e-3: continue
        tested += 1
        nearest_e = angle_of(Rm) <= min(angle_of(R[q].T @ Rm) for q in even if q != (0, 1, 2, 3)) + 1e-12
        if nearest_e != (l1 <= 1): bad += 1
    ok &= bad == 0 and tested > 2500
    vert = sorted(tuple(np.round(rodrigues(R[q]), 6)) for q in allp if cycle_type(q) == (4,))
    ok &= vert == sorted(tuple(np.round(np.array(v, float), 6)) for v in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)])
    print('ball: the identity cell is the octahedron |x|+|y|+|z|<=1 (mismatches %d of %d samples); its vertices are the six starred 4-cycles, the quarter-turns about x, y, z' % (bad, tested))
    # the loop: a third of a turn about the bond to hydrogen 1, +2pi/3, returns X0 with 2, 3, 4 cycled (an even relabelling)
    AX = np.array(X0[0]) / np.linalg.norm(X0[0])
    Rg = np.array(G.rot_axis(tuple(AX), 2 * math.pi / 3))
    Xg = [tuple(Rg @ np.array(x)) for x in X0]
    where = tuple(min(range(4), key=lambda j: np.linalg.norm(np.array(Xg[i]) - np.array(X0[j]))) for i in range(4))
    ok &= all(np.linalg.norm(np.array(Xg[i]) - np.array(X0[where[i]])) < 1e-9 for i in range(4)) and parity(list(where)) > 0 and where[0] == 0
    g = [q for q in even if np.allclose(R[q], Rg)]
    ok &= len(g) == 1
    cyclic = np.allclose(np.abs(Rg), np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])) or np.allclose(np.abs(Rg), np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]))
    ok &= cyclic
    # face gluing: the face x+y+z=1 (midplane to R_g) goes to x+y+z=-1 under R -> R R_g^-1, and the vertex R_x(pi/2) to
    # a quarter-turn about y or z (a third of a twist), not to R_x(-pi/2) (straight across)
    for _ in range(100):
        a, b = rng.random(2)
        if a + b > 1: a, b = 1 - a, 1 - b
        rod = np.array([a, b, 1 - a - b]); ang = 2 * math.atan(np.linalg.norm(rod))
        Rm = np.array(G.rot_axis(tuple(rod / np.linalg.norm(rod)), ang))
        ok &= abs(angle_of(Rm) - angle_of(Rg.T @ Rm)) < 1e-9                       # equidistant from e and R_g
        rod2 = rodrigues(Rm @ Rg.T)
        ok &= abs(rod2.sum() + 1) < 1e-9 and abs(angle_of(Rm @ Rg.T) - ang) < 1e-9  # lands on x+y+z=-1 at the same angle
    vx = np.array(G.rot_x(math.pi / 2)) @ Rg.T
    img = rodrigues(vx)
    ok &= abs(np.linalg.norm(img) - 1) < 1e-9 and abs(img[0]) < 1e-9
    print('loop: X0 returns with', [w + 1 for w in where], '| R_g cycles the axes x, y, z:', cyclic,
          '| the face x+y+z=1 is glued to x+y+z=-1 by R -> R R_g^-1; the vertex R_x(pi/2) lands at', np.round(img, 6).tolist())
    # the proper rotations T: 1E(h) = exp(i theta_h), theta_h the signed angle of R_h about the axis from C to the fixed H
    w = cmath.exp(2j * math.pi / 3)
    def oneE(q):
        if q == (0, 1, 2, 3): return 1
        if cycle_type(q) == (2, 2): return 1
        fixed = [i for i in range(4) if q[i] == i][0]
        axis = np.array(X0[fixed]) / np.linalg.norm(X0[fixed])
        Rh = R[q]
        s = (Rh[2, 1] - Rh[1, 2], Rh[0, 2] - Rh[2, 0], Rh[1, 0] - Rh[0, 1])   # 2 sin(theta) * axis
        return w if np.dot(s, axis) > 0 else w ** 2
    hom = all(abs(oneE(a) * oneE(b) - oneE(tuple(a[b[i]] for i in range(4)))) < 1e-9 for a in even for b in even)
    ok &= hom and abs(oneE(g[0]) - w) < 1e-9
    TT = {'A': lambda q: 1, '1E': oneE, '2E': lambda q: oneE(q).conjugate() if isinstance(oneE(q), complex) else oneE(q),
          'T': lambda q: np.trace(R[q])}
    spin_t = {k: round((sum(chi[q] * complex(f(q)).conjugate() for q in even) / 12).real) for k, f in TT.items()}
    ok &= spin_t == {'A': 5, '1E': 1, '2E': 1, 'T': 3}
    print('1E is a homomorphism of the even relabellings:', hom, '| 1E(g) = omega | spin under T:', spin_t)
    # packets: shares of a real packet with class overlaps c3 (third-turns), c2 (half-turns) against the closed forms
    for c3, c2 in ((0.0, 0.0), (0.2, 0.05), (0.37, 0.11)):
        cq = {q: 1.0 if q == (0, 1, 2, 3) else (c3 if cycle_type(q) == (3, 1) else c2) for q in even}
        share = {k: (complex(f((0, 1, 2, 3))).real / 12) * sum(complex(f(q)).conjugate() * cq[q] for q in even) for k, f in TT.items()}
        closed_form = {'A': (1 + 8 * c3 + 3 * c2) / 12, '1E': (1 - 4 * c3 + 3 * c2) / 12, '2E': (1 - 4 * c3 + 3 * c2) / 12, 'T': 3 * (1 - c2) / 4}
        ok &= all(abs(share[k] - closed_form[k]) < 1e-12 for k in TT) and abs(sum(share.values()) - 1) < 1e-12
    print('packet shares agree with the closed forms (1+8c3+3c2)/12, (1-4c3+3c2)/12, 3(1-c2)/4, summing to 1')
    # isomers: Gamma_rot x Gamma_nuc contains A (the even relabellings carry chi_stat = +1)
    iso = []
    for r, fr in TT.items():
        for nu, fn in TT.items():
            m = round((sum(complex(fr(q)) * complex(fn(q)) for q in even) / 12).real)
            if m: iso.append((r, nu, int(round(complex(fr((0, 1, 2, 3))).real)), spin_t[nu]))
    ok &= iso == [('A', 'A', 1, 5), ('1E', '2E', 1, 1), ('2E', '1E', 1, 1), ('T', 'T', 3, 3)]
    print('isomers (Gamma_rot, Gamma_nuc, d, m_nuc):', iso)
    # entanglement: invariants of T x T, the physical state of a T level with one T copy of the spins
    A = np.vstack([np.kron(R[q], R[q]) - np.eye(9) for q in even])
    _, s, vt = np.linalg.svd(A)
    null = vt[np.sum(s > 1e-9):]
    rank = np.linalg.matrix_rank(null[0].reshape(3, 3), tol=1e-9) if len(null) == 1 else None
    ok &= len(null) == 1 and rank == 3
    print('invariant vectors in T x T:', len(null), '| Schmidt rank of the physical state:', rank, '| entangled share 3 x 3 of 16 = 9/16')
    # monodromy: |T / ker Gamma|
    mono = {k: round(12 / sum(1 for q in even if abs(complex(f(q)) - complex(f((0, 1, 2, 3)))) < 1e-9)) for k, f in TT.items()}
    ok &= mono == {'A': 1, '1E': 3, '2E': 3, 'T': 12}
    print('monodromy group orders |T / ker Gamma|:', mono)
    # against the emitted numbers
    S = json.load(open(os.path.join(ROOT, 'figures', 'data', 'summary.json')))
    S0, S1 = S['methane_loop'], S['methane_rigid']
    ok &= S0['nucleus_i_sits_where_j_was'] == [w + 1 for w in where] and S0['axis'] == 'C-H1'
    ok &= S1['spin_Td'] == spin_td and S1['weights'] == weights and S1['spin_T'] == spin_t
    ok &= [tuple(x) for x in S1['isomers']] == iso and S1['entangled_fraction'] == [9, 16] and S1['monodromy_orders'] == mono
    ok &= S1['gamma_3N'] == gamma_3n and S1['gamma_vib'] == gamma_vib and S1['J_table'] == j_table
    print('summary.json agrees:', ok)
    if not ok:
        sys.exit('verify_methane: mismatch')
    print('verify_methane: all checks pass')

if __name__ == '__main__':
    main()
