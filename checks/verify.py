"""Checks of the group facts and of what the methylamine plates draw (python3 stdlib only).

Run from the repository root:  python3 checks/verify.py
Exact: the orders of H, G6, G12, B, S (2, 6, 12, 24, 240) and their cosets; [H,B] has 10 subgroups and 17 covers,
with three order-12 extensions of G6, any two generating B; btb = t^-1 and bt = t^2 b; the action of t and b on
G6/H and C[G6/H] = A1 + E; the proton-spin weights 12, 4, 8 (even parity) and 4, 12, 8 (odd) with chi_stat trivial
on G6; the KRb groups and the 16 subgroups of S; P_st = P_t P_s; C[G6] = C[G6/H] + A2 + E, with A2 + E induced from
the sign of H; on the (tau, iota) cylinder a shape is fixed by a conjugate of b on the six lines tau in (pi/3) Z, by a
conjugate of (23)* at the six points (pi/6 + k pi/3, 0), and by E alone elsewhere; the torsion species.
What the plates draw (figures/data/numbers.tex): figure 8's spin counts, figure 12's kappa rows and torsion species,
figure 5's version labels and figure 10's chart cells, all exactly; figure 9's widths (120 x degrees), the exponents
of c = exp(-1/(8 x^2)) under its bars, its shares against (1 + 2c)/3, its curve points, and the bars' split at the
no-overlap share 1/3 (dimension times copies of A1 in C[G6/H], over the three versions). summary.json: the intervals, the normal frame (M_tu = diag(1, -1)), the version positions and
the residual between X and bX. Geometry, in floating point: no element of G12 but E fixes X0.
verify_symbolic.py proves the family's actions and the Gaussian formulas exactly, verify.g repeats the group facts
in GAP. A failed assertion stops the build.
"""
import json, math, os, re, sys
from fractions import Fraction as Fr
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'compute'))
import geometry as G
import groups as Q

DATA = os.path.join(ROOT, 'figures', 'data')
NUM = dict(re.findall(r'\\def\\([A-Za-z]+)\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}', open(os.path.join(DATA, 'numbers.tex')).read()))
S = json.load(open(os.path.join(DATA, 'summary.json')))

def items(name):
    """The comma-separated items of a list macro, each split at the slashes outside braces."""
    out, depth, cur = [], 0, ''
    for ch in NUM[name] + ',':
        depth += (ch == '{') - (ch == '}')
        if ch == ',' and depth == 0: out.append(cur); cur = ''
        else: cur += ch
    fields = []
    for it in out:
        parts, depth, cur = [], 0, ''
        for ch in it + '/':
            depth += (ch == '{') - (ch == '}')
            if ch == '/' and depth == 0: parts.append(cur[1:-1] if cur.startswith('{') and cur.endswith('}') else cur); cur = ''
            else: cur += ch
        fields.append(parts)
    return fields

# ---- groups
assert [len(g) for g in (Q.H, Q.G6, Q.G12, Q.B, Q.S)] == [2, 6, 12, 24, 240]
assert Q.mul(Q.b, Q.mul(Q.t, Q.b)) == Q.inv(Q.t)                           # b t b = t^-1
assert Q.mul(Q.b, Q.t) == Q.mul(Q.mul(Q.t, Q.t), Q.b) != Q.mul(Q.t, Q.b)   # bt = t^2 b != tb
assert Q.close([Q.b, Q.t, Q.tu], Q.N) == Q.G12                              # G12 = <G6, tu>
assert len(Q.cosets(Q.G6, Q.H)) == 3 and len(Q.cosets(Q.G12, Q.H)) == 6
assert len(Q.S) // len(Q.G12) == 20 and len(Q.S) // len(Q.H) == 120
subs, names, cov = Q.interval_data()
assert len(subs) == 10 and len(cov) == 17
assert sorted(len(K) for K in subs) == [2, 4, 4, 4, 6, 8, 12, 12, 12, 24]
order12 = [K for K in subs if len(K) == 12 and Q.G6 <= K]
assert len(order12) == 3 and all(Q.close(list(a | b), Q.N) == Q.B for i, a in enumerate(order12) for b in order12[i+1:])
_, act_t, act_b = Q.coset_action_table()
assert act_t == [1, 2, 0] and act_b == [0, 2, 1]        # t cycles the three cosets, b fixes H and swaps tH and t^2H
# C[G6/H]: the permutation character counts the fixed cosets
size = {'E': 1, 't': 2, 'b': 3}
perm_char = {'E': 3, 't': sum(1 for i, j in enumerate(act_t) if i == j), 'b': sum(1 for i, j in enumerate(act_b) if i == j)}
local = {n: Fr(sum(size[c] * chi[c] * perm_char[c] for c in size), 6) for n, chi in Q.G6_CHARS.items()}
assert local == {'A1': 1, 'A2': 0, 'E': 1}
w = Q.spin_weights()
assert w == {1: {'A1': 12, 'A2': 4, 'E': 8}, -1: {'A1': 4, 'A2': 12, 'E': 8}}
assert all(v['A1'] + v['A2'] + 2 * v['E'] == 32 for v in w.values())
assert all(Q.sign(g[0][:5]) == 1 for g in Q.G6)        # chi_stat is trivial on G6
Sk, Gin, GK, GRb, _, _, _ = Q.krb_groups()
assert (len(Sk), len(Gin), len(GK), len(GRb)) == (8, 4, 4, 4)
assert all(len(Q.close(list(a | b), 4)) == 8 for a, b in [(Gin, GK), (Gin, GRb), (GK, GRb)])
assert len(Q.interval(Q.close([], 4), Sk, 4)) == 16
def pm(g):                                              # P_s with the inversion sign
    return [[(-1 if g[1] else 1) * v for v in row] for row in G.perm_matrix(g[0], Q.N)]
def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
assert all(pm(Q.mul(s, t_)) == mm(pm(t_), pm(s)) for s in (Q.t, Q.u, Q.b) for t_ in (Q.t, Q.u, Q.b, Q.tu))
print('groups, cosets, intervals, C[G6/H] = A1 + E, spin weights, KRb, P_st = P_t P_s: ok')

# ---- the base point: induced spaces and the point groups of shapes (Sections 7 and 10)
def induced(K, chi_K):
    """Multiplicities in G6 of the representation induced from the character chi_K of K <= G6."""
    conj = lambda x, g: Q.mul(Q.inv(x), Q.mul(g, x))
    ind = {g: Fr(sum(chi_K[conj(x, g)] for x in Q.G6 if conj(x, g) in K), len(K)) for g in Q.G6}
    return {n: Fr(sum(ch[Q.g6_class(g)] * ind[g] for g in Q.G6), 6) for n, ch in Q.G6_CHARS.items()}
E7 = Q.E(Q.N)
assert induced({E7}, {E7: 1}) == {'A1': 1, 'A2': 1, 'E': 2}                 # C[G6]
assert induced(Q.H, {E7: 1, Q.b: 1}) == {'A1': 1, 'A2': 0, 'E': 1}         # C[G6/H]
assert induced(Q.H, {E7: 1, Q.b: -1}) == {'A1': 0, 'A2': 1, 'E': 1}        # induced from the sign of H
# a shape (tau, iota) moves as tau -> e tau + k pi/3, iota -> (-1)^k iota, with t, u, b as verify_symbolic.py proves
gen = {Q.t: (1, 2), Q.u: (1, 3), Q.b: (-1, 0)}
shape = {E7: (1, 0)}; todo = [E7]
while todo:
    g = todo.pop()
    for h, (e2, k2) in gen.items():
        e1, k1 = shape[g]; gh = Q.mul(g, h); lab = (e1 * e2, (e1 * k2 + k1) % 6)   # A_gh = A_g o A_h
        if gh not in shape:
            shape[gh] = lab; todo.append(gh)
        assert shape[gh] == lab                                                # well defined on the group
assert set(shape) == Q.G12 and len(set(shape.values())) == 12                  # a faithful action
def stabilizer(tau, iota):                                                     # tau in units of pi, iota in units of iota0
    return frozenset(g for g, (e, k) in shape.items() if (e * tau + Fr(k, 3) - tau) % 2 == 0 and (-1) ** k * iota == iota)
cls = lambda g: frozenset(Q.mul(x, Q.mul(g, Q.inv(x))) for x in Q.G12)
c23 = Q.cyc(Q.N, (2, 3), star=1)
assert c23 in Q.G12 and c23 not in cls(Q.b) and len(cls(Q.b)) == len(cls(c23)) == 3
for tau in [Fr(j, 12) for j in range(24)] + [Fr(1, 7), Fr(9, 7)]:
    for iota in (Fr(-1), Fr(-1, 2), Fr(0), Fr(1, 3), Fr(1)):
        K = stabilizer(tau, iota) - {E7}
        if tau % Fr(1, 3) == 0:
            assert len(K) == 1 and next(iter(K)) in cls(Q.b), (tau, iota)        # the six lines tau in (pi/3) Z
        elif iota == 0 and (tau - Fr(1, 6)) % Fr(1, 3) == 0:
            assert len(K) == 1 and next(iter(K)) in cls(c23), (tau, iota)        # the six points (pi/6 + k pi/3, 0)
        else:
            assert not K, (tau, iota)
print('C[G6] = C[G6/H] + (A2 + E), the sign of H induces A2 + E; shapes fixed by conjugates of b, of (23)*, or by E alone: ok')

# ---- the torsion species: (cos m tau, sin m tau) carries t with trace 2 cos(2 pi m/3) (2 when 3 | m, else -1), b with 0
species = [['0', 'A1', 'none']]
for m in range(1, 7):
    chi = {'E': 2, 't': 2 if m % 3 == 0 else -1, 'b': 0}
    mult = {n: Fr(sum(size[c] * ch[c] * chi[c] for c in size), 6) for n, ch in Q.G6_CHARS.items()}
    assert mult in ({'A1': 0, 'A2': 0, 'E': 1}, {'A1': 1, 'A2': 1, 'E': 0})
    species.append([str(m)] + (['E', 'E'] if mult['E'] else ['A1', 'A2']))   # cos is even under b, sin odd
assert items('torsionSpecies') == species
assert [[str(m), c, s or 'none'] for m, c, s in S['torsion_species']] == species

# ---- what the plates draw
assert [int(NUM[k]) for k in ('spinAone', 'spinAtwo', 'spinE', 'spinDim')] == [12, 4, 8, 32]
assert Fr(NUM['rhoIll']) == Fr(1, 2)
assert [[r[0], [Fr(x) for x in r[1].split(',')]] for r in items('kappaRows')] == \
       [[str(K), [m + Fr(K, 2) for m in range(-3, 4)]] for K in (-2, -1, 0, 1, 2)]
words = {'E': [], 't': ['t'], 't^2': ['t', 't'], 'u': ['u'], 'tu': ['t', 'u'], 't^2u': ['t', 't', 'u']}
for row in items('versionList'):
    g = Q.E(Q.N)
    for x in words[row[1]]:
        g = Q.mul(g, {'t': Q.t, 'u': Q.u}[x])
    assert [int(v) for v in row[2:]] == [g[0][j] + 1 for j in range(7)] and row[0] == {'E': 'H'}.get(row[1], row[1] + 'H')
# the chart cells: U = [-1/3, 1/3] x [0, 1] (tau in units of pi, iota in units of iota0) moved by each coset
# representative, t: tau + 2/3, u: tau + 1 and iota -> -iota, reduced to [-1, 1] and split at the seam
cells = []
for word in words.values():
    shift, side = Fr(0), 1
    for x in word:
        shift += Fr(2, 3) if x == 't' else 1; side *= 1 if x == 't' else -1
    lo = (shift - Fr(1, 3) + 1) % 2 - 1; hi = lo + Fr(2, 3)
    e = (0, 1) if side > 0 else (-1, 0)
    cells += [(lo, min(hi, Fr(1))) + e] + ([(Fr(-1), hi - 2) + e] if hi > 1 else [])
emitted = sorted((Fr(a), Fr(b), int(c), int(d)) for a, b, c, d in items('chartCells'))
assert len(emitted) == len(cells) and all(abs(x - y) < 1e-4 for p, q in zip(emitted, sorted(cells)) for x, y in zip(p, q))
# figure 9: each row's width is 120 x degrees (d = 120), its c = e^{-1/(8 x^2)} has the exponent drawn under its bar, and
# its share is (1 + 2c)/3. The bars split at the no-overlap share 1/3, one site of three, and the curve points are the rows.
rows9 = items('packetRows')
for x, deg, share, en, ed in rows9:
    xf = Fr(x)
    assert Fr(deg) == 120 * xf and Fr(int(en), int(ed)) == Fr(1, 8) / xf ** 2
    assert abs(float(share) - (1 + 2 * math.exp(-float(Fr(int(en), int(ed))))) / 3) < 5e-5
assert local['A1'] == 1 and abs(float(NUM['packetFloor']) - float(Fr(1, 1) * local['A1'] / len(Q.cosets(Q.G6, Q.H)))) < 5e-5
marks = re.findall(r'\(([\d.]+),([\d.]+)\)', NUM['packetMarks'])
assert [(round(float(Fr(x)), 4), share) for x, _, share, _, _ in rows9] == [(float(a), b) for a, b in marks]
with open(os.path.join(DATA, 'gaussian-shares.dat')) as f:
    for x, a, e in (l.split() for l in f.read().strip().splitlines()[1:]):
        c = math.exp(-1 / (8 * float(x) ** 2)) if float(x) > 0 else 0.0   # the x = 0 row carries the limit c -> 0
        assert abs(float(a) - (1 + 2 * c) / 3) < 1e-6 and abs(float(e) - 2 * (1 - c) / 3) < 1e-6
print('numbers.tex: spin counts, kappa rows, torsion species, version labels, chart cells, packet widths, exponents, shares, marks and the 1/3 split: ok')

# ---- summary.json and the geometry of the reference
assert S['interval'] == {'subgroups': 10, 'covers': 17, 'versions': [1, 2, 2, 2, 3, 4, 6, 6, 6, 12]}
assert S['interval_HS'] == {'subgroups': 36, 'covers': 73, 'in_bond_interval': 10}
assert S['spin_weights'] == {str(k): v for k, v in w.items()}
M = S['normal_frame']['M']
assert abs(M[0][0] - 1) < 1e-9 and abs(M[1][1] + 1) < 1e-9 and abs(M[0][1]) < 1e-9 and abs(M[1][0]) < 1e-9
assert S['normal_frame']['resid'] < 1e-9 and S['normal_frame']['ortho'] < 1e-12 and S['normal_frame']['tangent_dim'] == 5
vpos = {c: (ta, io) for c, g, ta, io in S['versions']}
assert abs(vpos['tH'][0] - 2/3) < 1e-6 and abs(vpos['tuH'][0] - 5/3) < 1e-6 and vpos['tuH'][1] < 0 < vpos['H'][1]
assert S['generic_X_bX_best_rotation_residual'] > 0.1          # X and bX are not related by a rotation
X0 = G.methylamine(0.0, G.MLA_FRAME['iota0'])
dmin = min(G.config_distance(G.apply_perm_inversion(X0, g[0], g[1]), X0) for g in Q.G12 if g != Q.E(Q.N))
assert dmin > 1.0
print(f'summary.json: intervals, normal frame, versions, X and bX; G12 acts freely near X0 (nearest image {dmin:.2f} A): ok')
print('verify.py: all checks passed')
