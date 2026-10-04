"""Verification tests for the figure data (python3 stdlib only).

Run from the repository root:  python3 checks/verify.py
Exact (integer / rational) checks come first; the geometric model checks use
floating point and report their residuals.  Any assertion failure is a build
failure.  checks/verify_symbolic.py repeats the trigonometric identities and
the Gaussian formulas symbolically (needs sympy); checks/verify.g repeats the
group facts in GAP.
"""
import json, math, os, sys
from fractions import Fraction as Fr
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'compute'))
import geometry as G
import groups as Q

report = {}

# ------------------------------------------------------------ exact: groups --
assert [len(g) for g in (Q.H, Q.G6, Q.G12, Q.B, Q.S)] == [2, 6, 12, 24, 240]
assert Q.mul(Q.b, Q.mul(Q.t, Q.b)) == Q.inv(Q.t)                      # b t b = t^-1
assert Q.mul(Q.b, Q.t) == Q.mul(Q.mul(Q.t, Q.t), Q.b) != Q.mul(Q.t, Q.b)   # bt = t^2 b != tb
assert Q.close([Q.b, Q.t, Q.tu], Q.N) == Q.G12                         # G12 = <G6, tu>
assert len(Q.cosets(Q.G6, Q.H)) == 3 and len(Q.cosets(Q.G12, Q.H)) == 6
assert len(Q.S) // len(Q.G12) == 20 and len(Q.S) // len(Q.H) == 120
subs, names, cov = Q.interval_data()
assert len(subs) == 10 and len(cov) == 17
assert sorted(len(K) for K in subs) == [2, 4, 4, 4, 6, 8, 12, 12, 12, 24]
order12 = [K for K in subs if len(K) == 12 and Q.G6 <= K]
assert len(order12) == 3
assert all(Q.close(list(a | b), Q.N) == Q.B for i, a in enumerate(order12) for b in order12[i+1:])
names_, act_t, act_b = Q.coset_action_table()
assert act_t == [1, 2, 0] and act_b == [0, 2, 1]        # t cycles, b fixes H and swaps tH, t^2H
w = Q.spin_weights()
assert w == {1: {'A1': 12, 'A2': 4, 'E': 8}, -1: {'A1': 4, 'A2': 12, 'E': 8}}
for v in w.values():
    assert v['A1'] + v['A2'] + 2 * v['E'] == 32
# statistics character is trivial on G6 (all underlying permutations are even)
assert all(Q.sign(g[0][:5]) == 1 for g in Q.G6)
# KRb
Sk, Gin, GK, GRb, p12, p34, Es = Q.krb_groups()
assert (len(Sk), len(Gin), len(GK), len(GRb)) == (8, 4, 4, 4)
assert all(len(Q.close(list(a | b), 4)) == 8 for a, b in [(Gin, GK), (Gin, GRb), (GK, GRb)])
subsK = Q.interval(Q.close([], 4) | {Q.E(4)}, Sk, 4)
assert len(subsK) == 16
report['groups'] = 'orders, cosets, interval [H,B] (10 subgroups, 17 covers), relations, weights, KRb: ok'

# ---------------------------------------------- exact: permutation matrices --
def pm(g):   # column-permutation matrix as nested lists of ints with inversion sign
    P = G.perm_matrix(g[0], Q.N)
    return [[(-1 if g[1] else 1) * v for v in row] for row in P]
def mm(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
for s in (Q.t, Q.u, Q.b):
    for tt in (Q.t, Q.u, Q.b, Q.tu):
        assert pm(Q.mul(s, tt)) == mm(pm(tt), pm(s))                  # P_{st} = P_t P_s
report['perm_matrices'] = 'P_{st} = P_t P_s on generators: ok'

# ------------------------------------------- exact: coefficient basis (Q(sqrt)) --
# Work with squared norms and inner products as rationals: v1=(1,1,1)/s3, v2=(2,-1,-1)/s6, v3=(0,1,-1)/s2.
vecs = {'v1': ([1, 1, 1], 3), 'v2': ([2, -1, -1], 6), 'v3': ([0, 1, -1], 2)}
for a, (va, na) in vecs.items():
    for b_, (vb, nb) in vecs.items():
        ip = sum(x*y for x, y in zip(va, vb))        # inner product times sqrt(na nb)
        assert (ip == 0) if a != b_ else (Fr(ip, na) == 1)
# site action of t in the ordered basis (H, tH, t^2H): e1->e2->e3->e1 ; b: e2<->e3
T = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]; Bm = [[1, 0, 0], [0, 0, 1], [0, 1, 0]]
def mv(M, v): return [sum(M[i][j]*v[j] for j in range(3)) for i in range(3)]
assert mv(T, [1, 1, 1]) == [1, 1, 1] and mv(Bm, [1, 1, 1]) == [1, 1, 1]
# T v2 = -1/2 v2 + (sqrt3/2) v3  <=>  2 T v2 = -v2 + sqrt3 v3; compare coordinates with sqrt3 v3 = (0, sqrt3, -sqrt3)/sqrt2 ... use the unnormalised vectors:
# in unnormalised form: T(2,-1,-1) = (-1,2,-1) = -1/2 (2,-1,-1) + 3/2 (0,1,-1)
assert mv(T, [2, -1, -1]) == [-1, 2, -1]
assert all(Fr(x) == Fr(-1, 2)*a + Fr(3, 2)*c for x, a, c in zip([-1, 2, -1], [2, -1, -1], [0, 1, -1]))
# T(0,1,-1) = (-1,0,1) = -1/2 (2,-1,-1) - 1/2 (0,1,-1)
assert mv(T, [0, 1, -1]) == [-1, 0, 1]
assert all(Fr(x) == Fr(-1, 2)*a + Fr(-1, 2)*c for x, a, c in zip([-1, 0, 1], [2, -1, -1], [0, 1, -1]))
# with normalisation this is the rotation by 120 degrees: entries (-1/2, -sqrt3/2; sqrt3/2, -1/2)
assert mv(Bm, [2, -1, -1]) == [2, -1, -1] and mv(Bm, [0, 1, -1]) == [0, -1, 1]    # b = diag(1,-1) on E
report['coefficient_basis'] = 'orthonormal; t -> rotation by 2pi/3, b -> diag(1,-1) on the E pair: ok'

# -------------------------------------------------- exact: torsion species --
# U_t cos(m tau) = cos(m tau - 2 pi m/3) etc.  Character of t on {cos,sin} is 2cos(2 pi m/3) in {2,-1}.
for m in range(0, 7):
    if m == 0:
        continue
    chi_t = 2 if m % 3 == 0 else -1
    mults = {n: Fr(Q.G6_CHARS[n]['E']*2 + 2*Q.G6_CHARS[n]['t']*chi_t + 0, 6) for n in Q.G6_CHARS}
    if m % 3:
        assert mults == {'A1': 0, 'A2': 0, 'E': 1}
    else:
        assert mults == {'A1': 1, 'A2': 1, 'E': 0}
report['torsion_species'] = 'm=0: A1; m not divisible by 3: E pair; nonzero multiple of 3: A1 (cos) + A2 (sin): ok'

# ---------------------------------------------------- exact: Gaussian shares --
# sum rule and limits as rational identities in c
for c in (Fr(0), Fr(1, 3), Fr(1, 2), Fr(1)):
    assert Fr(1 + 2*c, 3) + Fr(2*(1 - c), 3) == 1
assert (Fr(1, 3), Fr(2, 3)) == (Fr(1 + 0, 3), Fr(2*(1 - 0), 3)) and (Fr(3, 3), Fr(0)) == (Fr(1 + 2, 3), Fr(0))
# numerical confirmation of c = exp(-d^2/(8 Delta^2)) for density variance Delta^2 (2D grid quadrature)
def overlap_numeric(d, Delta, n=400, L=8.0):
    h = 2*L/n; s = 0.0; nrm = 0.0
    for i in range(n):
        x = -L + (i + .5)*h
        for j in range(n):
            y = -L + (j + .5)*h
            g0 = math.exp(-(x*x + y*y)/(4*Delta*Delta))
            g1 = math.exp(-((x - d)**2 + y*y)/(4*Delta*Delta))
            s += g0*g1; nrm += g0*g0
    return s/nrm
for d, sg in ((1.0, 0.5), (1.0, 1.0), (2.0, 0.8)):
    assert abs(overlap_numeric(d, sg) - math.exp(-d*d/(8*sg*sg))) < 1e-6
report['gaussian'] = 'sum rule and limits exact; overlap formula confirmed by quadrature to 1e-6'

# ------------------------------------------------------- numeric: geometry --
X0 = G.methylamine(0.0, G.MLA_FRAME['eta0'])
assert max(abs(v) for v in G.mass_moment(X0, G.MLA_MASSES)) < 1e-12
bX0 = G.apply_perm_inversion(X0, Q.b[0], Q.b[1])
res_b = G.config_distance(bX0, G.rotate(G.rot_y(math.pi), X0))
assert res_b < 1e-12
# the (tau, eta) actions derived in docs/verification.md
eta0 = G.MLA_FRAME['eta0']
checks = {
    't': (Q.t, lambda a, e: (a + 2*math.pi/3, e), G.identity(3)),
    'u': (Q.u, lambda a, e: (a + math.pi, -e), G.rot_z(math.pi)),
    'b': (Q.b, lambda a, e: (-a, e), G.rot_y(math.pi)),
    'tu': (Q.tu, lambda a, e: (a - math.pi/3, -e), G.rot_z(math.pi)),
}
worst = 0.0
for name, (g, f, A) in checks.items():
    for a, e in ((0.0, eta0), (0.8, -0.3), (-2.4, 0.1), (3.0, eta0)):
        X = G.methylamine(a, e)
        gX = G.apply_perm_inversion(X, g[0], g[1])
        tp, ep = f(a, e)
        worst = max(worst, G.config_distance(gX, G.rotate(A, G.methylamine(tp, ep))))
assert worst < 1e-12
# the water configurations are centred and the rotation/relabelling actions commute
wX = G.water(0.2, 112.0)
assert max(abs(v) for v in G.mass_moment(wX, G.WATER_MASSES)) < 1e-12
Rw = G.rot_z(0.7)
sig = [1, 0, 2]
assert G.config_distance(G.rotate(Rw, G.apply_perm_inversion(wX, sig, 0)),
                         G.apply_perm_inversion(G.rotate(Rw, wX), sig, 0)) < 1e-12
# free action near the reference: no nontrivial element of G12 fixes X0
dmin = min(G.config_distance(G.apply_perm_inversion(X0, g[0], g[1]), X0) for g in Q.G12 if g != Q.E(Q.N))
assert dmin > 1.0
report['geometry'] = f'centring < 1e-12; b = R_y(pi) residual {res_b:.1e}; family actions residual {worst:.1e}; min |gX0-X0| = {dmin:.3f} A'

# ----------------------------------------- numeric: emitted data consistency --
with open(os.path.join(ROOT, 'figures', 'data', 'summary.json')) as f:
    S = json.load(f)
assert S['spin_weights'] == {'1': {'A1': 12, 'A2': 4, 'E': 8}, '-1': {'A1': 4, 'A2': 12, 'E': 8}}
assert S['local_copies'] == {'A1': 1, 'A2': 0, 'E': 1}
assert S['interval'] == {'subgroups': 10, 'covers': 17, 'versions': [1, 2, 2, 2, 3, 4, 6, 6, 6, 12]}
assert S['interval_HS']['subgroups'] == 36 and S['interval_HS']['covers'] == 73 and S['interval_HS']['in_bond_interval'] == 10
M = S['normal_frame']['M']
assert abs(M[0][0] - 1) < 1e-9 and abs(M[1][1] + 1) < 1e-9 and abs(M[0][1]) < 1e-9 and abs(M[1][0]) < 1e-9
assert S['normal_frame']['resid'] < 1e-9 and S['normal_frame']['ortho'] < 1e-12
assert S['normal_frame']['tangent_dim'] == 5
vpos = {c: (ta, et) for c, g, ta, et in S['versions']}
assert abs(vpos['tH'][0] - 2/3) < 1e-6 and abs(vpos['tuH'][0] - 5/3) < 1e-6 and vpos['tuH'][1] < 0 < vpos['H'][1]
assert S['generic_X_bX_best_rotation_residual'] > 0.1          # X and bX are not rotation related
assert [m for m, _, _ in S['torsion_species']] == list(range(7))
with open(os.path.join(ROOT, 'figures', 'data', 'gaussian-shares.dat')) as f:
    rows = [l.split() for l in f.read().strip().splitlines()[1:]]
for x, a, e, z in rows:
    c = math.exp(-1/(8*float(x)**2)) if float(x) > 0 else 0.0   # the x = 0 row carries the limit c -> 0
    assert abs(float(a) - (1 + 2*c)/3) < 1e-6 and abs(float(e) - 2*(1 - c)/3) < 1e-6 and float(z) == 0
report['emitted_data'] = 'summary.json and gaussian-shares.dat agree with recomputation'

print(json.dumps(report, indent=1))
print('verify.py: all checks passed')
