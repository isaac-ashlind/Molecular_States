"""The closing pages: rigid methane run through the pipeline and set against Albert, Kubischta, Lemeshko and Liu
(arXiv:2403.04572v4), in exact arithmetic (integers, rationals and the cyclotomic field Q(omega), omega = e^{2 pi i/3}),
and compared with what compute/make_data.py emits (summary.json, keys methane_rigid and methane_loop).

  geometry   X0 is the regular tetrahedron with its hydrogens on the body diagonals d1..d4 of a cube (the C2 frame).
             Each element of Td(M), an even relabeling or an odd one starred, carries X0 into its orbit by a proper
             rotation R_h, a signed permutation matrix, and no other relabeling does. The 24 close as the cube group O
             (1, 6 quarter-turns, 8 third-turns, 9 half-turns), the 12 even as T, and R_{h1h2} = R_{h2} R_{h1}.
             |S| = 48 and |S/H| = 2.
  Table A    per class the turn of R_h and its axis, chi_stat, chi_-, chi_spin and chi_3N
  spin       (C^2)^{x4} x C^1 under Td(M) is 5 A1 + E + 3 T2, under T 5 A + 1E + 2E + 3 T; the physical states by parity
  vibrations the fifteen displacements are A1 + E + T1 + 3 T2, and without translations and rotations the nine normal
             ones are A1 + E + 2 T2
  J ladder   D^J restricted to H for J = 0..6, from the exact turning angles, and the physical states of each J by parity
  T          C[T] = A + 1E + 2E + 3T; the species of Td(M) on T; 1E is a homomorphism with 1E(g) = omega; the ring K
             under g carries A when 3 | K, 1E when K = 1 and 2E when K = 2 mod 3
  isomers    Gamma_rot x Gamma_nuc contains A for (A, A), (1E, 2E), (2E, 1E), (T, T), the same when 1E and 2E exchange
  entangled  the invariant of T x T has Schmidt rank 3, so 9 of the 16 spin states
  monodromy  |T / ker Gamma| = 1, 3, 3, 12; the third-turn about the bond to hydrogen 1 returns X0 with 2, 3, 4 cycled
  ball       in Rodrigues coordinates the cell of the identity among the twelve rotations of T is the octahedron
             |x| + |y| + |z| <= 1, its vertices the quarter-turns and the face centers of the cube of third-turns; the
             face x + y + z = 1 is glued to x + y + z = -1 by R -> R R_g^-1 with a third of a twist
  packets    the shares (1 + 8 c3 + 3 c2)/12, (1 - 4 c3 + 3 c2)/12, 3 (1 - c2)/4 as exact polynomials in c3 and c2
Exit status nonzero on any mismatch.
"""
import itertools, json, math, os, sys
from fractions import Fraction as Fr
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'compute'))
import geometry as G

class W:
    """a + b omega in Q(omega), omega = e^{2 pi i/3}, omega^2 = -1 - omega."""
    def __init__(self, a, b=0):
        if isinstance(a, W): a, b = a.a, a.b
        self.a, self.b = Fr(a), Fr(b)
    def __add__(self, o): o = o if isinstance(o, W) else W(o); return W(self.a + o.a, self.b + o.b)
    __radd__ = __add__
    def __mul__(self, o):
        o = o if isinstance(o, W) else W(o)
        return W(self.a * o.a - self.b * o.b, self.a * o.b + self.b * o.a - self.b * o.b)
    __rmul__ = __mul__
    def conj(self): return W(self.a - self.b, -self.b)
    def __eq__(self, o): o = o if isinstance(o, W) else W(o); return self.a == o.a and self.b == o.b
    def rational(self): assert self.b == 0; return self.a
OMEGA = W(0, 1)

# ---- small exact linear algebra on lists of Fractions or ints
def mat(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def tr(A): return sum(A[i][i] for i in range(len(A)))
def T_(A): return [list(r) for r in zip(*A)]
def vec(A): return (A[2][1] - A[1][2], A[0][2] - A[2][0], A[1][0] - A[0][1])
def rref(A):
    """Row-reduced echelon form over Q; returns (rows, pivot columns)."""
    A = [[Fr(x) for x in r] for r in A]; piv = []; r = 0
    for c in range(len(A[0])):
        k = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if k is None: continue
        A[r], A[k] = A[k], A[r]; A[r] = [x / A[r][c] for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0: A[i] = [x - A[i][c] * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    return A[:r], piv
def rank(A): return len(rref(A)[1])
def inverse(A):
    R, piv = rref([list(r) + [int(i == k) for k in range(len(A))] for i, r in enumerate(A)])
    assert piv == list(range(len(A)))
    return [row[len(A):] for row in R]
def det3(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
def nullspace(A):
    R, piv = rref(A); n = len(A[0]); out = []
    for f in (c for c in range(n) if c not in piv):
        v = [Fr(0)] * n; v[f] = Fr(1)
        for row, p in zip(R, piv): v[p] = -row[f]
        out.append(v)
    return out
def solve(B, y):
    """Coordinates c with B c = y for a full-column-rank B (columns are basis vectors)."""
    R, piv = rref([list(r) + [v] for r, v in zip(B, y)])
    assert len(B[0]) not in piv and piv == list(range(len(B[0])))
    return [row[-1] for row in R]

def parity(p): return (-1) ** sum(1 for i in range(4) for j in range(i + 1, 4) if p[i] > p[j])
def cycles(p, n=4):
    seen, lens = set(), []
    for i in range(n):
        if i in seen: continue
        k, j = 0, i
        while j not in seen: seen.add(j); j = p[j]; k += 1
        lens.append(k)
    return tuple(sorted(lens, reverse=True))

D = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]     # the bond directions, hydrogen i on d_i

def plain(x):
    """Fractions and nested dicts and tuples for printing: integers as integers, the rest as p/q."""
    if isinstance(x, dict): return {k: plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return type(x)(plain(v) for v in x)
    if isinstance(x, Fr): return int(x) if x.denominator == 1 else f'{x.numerator}/{x.denominator}'
    return x

def main():
    ok = True
    # ---- geometry: the emitted X0 is this tetrahedron, the carbon at the center
    X0 = G.methane_c2(); s = math.dist(X0[0], X0[4]) / math.sqrt(3)
    ok &= all(abs(X0[i][k] - s * D[i][k]) < 1e-12 for i in range(4) for k in range(3)) and max(map(abs, X0[4])) < 1e-12
    # h.X0 = X0 P_h: column i receives eps x_{sigma^-1(i)}, eps = -1 when starred; R_h solves R d_i = eps d_{sigma^-1(i)}
    D3inv = inverse(T_([list(d) for d in D[:3]]))
    R, proper = {}, {}
    for p in itertools.permutations(range(4)):
        for star in (False, True):
            eps = -1 if star else 1
            img = [tuple(eps * c for c in D[p.index(i)]) for i in range(4)]
            Rm = mat(T_([list(v) for v in img[:3]]), D3inv)
            is_rot = mat(Rm, [[x] for x in D[3]]) == [[x] for x in img[3]] and mat(Rm, T_(Rm)) == I3 and det3(Rm) == 1
            proper[(p, star)] = is_rot
            if is_rot: R[p] = [[int(x) for x in r] for r in Rm]
    ok &= all(proper[(p, star)] == (star == (parity(p) < 0)) for p, star in proper)   # exactly Td(M)
    allp = sorted(R); even = [p for p in allp if parity(p) > 0]; e = (0, 1, 2, 3)
    cosang = {p: Fr(tr(R[p]) - 1, 2) for p in allp}                                     # cos of the turn of R_h
    turn = {Fr(1): 0, Fr(0): 90, Fr(-1, 2): 120, Fr(-1): 180}
    ok &= sorted(turn[cosang[p]] for p in allp) == [0] + [90] * 6 + [120] * 8 + [180] * 9
    comp = lambda a, b: tuple(a[b[i]] for i in range(4))                               # (ab)(i) = a(b(i))
    ok &= all(comp(a, b) in R for a in allp for b in allp) and all(comp(a, b) in even for a in even for b in even) and len(even) == 12
    rev = all(R[comp(a, b)] == mat(R[b], R[a]) for a in allp for b in allp)
    fwd = all(R[comp(a, b)] == mat(R[a], R[b]) for a in allp for b in allp)
    ok &= rev and not fwd and 2 * math.factorial(4) == 48 and len(allp) == 24
    print('geometry: Td(M) carries X0 into its orbit by 24 signed permutation matrices (O), the 12 even ones T;',
          'R_{h1h2} = R_{h2} R_{h1}:', rev, '| |S| = 48, |H| = 24, |S/H| = 2')

    # ---- class functions on Td(M): Table A
    def axis(p):                                           # a half-turn's axis from R + 1 = 2 n n^T, any other's from R - R^T
        if cosang[p] != -1: return vec(R[p])
        return next(c for c in (tuple(R[p][i][k] + (i == k) for i in range(3)) for k in range(3)) if any(c))
    def axis_kind(p):
        if p == e: return ''
        nz = sorted(abs(x) for x in axis(p) if x)
        return {1: 'x, y, z', 2: 'cube edge', 3: 'bond'}[len(nz)] if len(set(nz)) == 1 else '?'
    def spin_trace(p):                                     # fixed basis states of the relabeling on (C^2)^{x4} x C^1
        return sum(1 for idx in range(16) if all(((idx >> p[i]) & 1) == ((idx >> i) & 1) for i in range(4)))
    def action15(p):                                       # dX -> R_h^T dX P_h on the 3 x 5 displacements, P_h = eps P_sigma
        eps = parity(p); P = [[eps if (i < 4 and p[i] == j) or (i == j == 4) else 0 for j in range(5)] for i in range(5)]
        cols = []
        for i in range(5):
            for a in range(3):
                E = [[1 if (r == a and c == i) else 0 for c in range(5)] for r in range(3)]
                M = mat(mat(T_(R[p]), E), P); cols.append([M[r][c] for c in range(5) for r in range(3)])
        return T_(cols)
    chi_spin = {p: spin_trace(p) for p in allp}
    A15 = {p: action15(p) for p in allp}
    chi_3n = {p: tr(A15[p]) for p in allp}
    table_a = {}
    for p in allp:
        table_a.setdefault(cycles(p), set()).add((parity(p) < 0, turn[cosang[p]], axis_kind(p), parity(p), parity(p), chi_spin[p], chi_3n[p]))
    ok &= table_a == {(1, 1, 1, 1): {(False, 0, '', 1, 1, 16, 15)}, (3, 1): {(False, 120, 'bond', 1, 1, 4, 0)},
                      (2, 2): {(False, 180, 'x, y, z', 1, 1, 4, -1)}, (4,): {(True, 90, 'x, y, z', -1, -1, 2, -1)},
                      (2, 1, 1): {(True, 180, 'cube edge', -1, -1, 8, 3)}}
    print('Table A (star, turn, axis, chi_stat, chi_-, chi_spin, chi_3N):', table_a)

    # ---- species of Td(M) (by class), the orthogonality relations, the spin decomposition and the pairing by parity
    TD = {'A1': {(1,1,1,1): 1, (3,1): 1, (2,2): 1, (4,): 1, (2,1,1): 1},
          'A2': {(1,1,1,1): 1, (3,1): 1, (2,2): 1, (4,): -1, (2,1,1): -1},
          'E':  {(1,1,1,1): 2, (3,1): -1, (2,2): 2, (4,): 0, (2,1,1): 0},
          'T1': {(1,1,1,1): 3, (3,1): 0, (2,2): -1, (4,): 1, (2,1,1): -1},
          'T2': {(1,1,1,1): 3, (3,1): 0, (2,2): -1, (4,): -1, (2,1,1): 1}}
    def td(f): return {k: Fr(sum(f(p) * v[cycles(p)] for p in allp), 24) for k, v in TD.items()}
    ok &= all(td(lambda p, v=v: v[cycles(p)])[k] == (k == g) for g, v in TD.items() for k in TD)
    spin_td = td(lambda p: chi_spin[p])
    ok &= spin_td == {'A1': 5, 'A2': 0, 'E': 1, 'T1': 0, 'T2': 3}
    weights = {tag: td(lambda p, par=par: chi_spin[p] * parity(p) * (par if parity(p) < 0 else 1)) for par, tag in ((1, 'Even'), (-1, 'Odd'))}
    ok &= weights == {'Even': {'A1': 0, 'A2': 5, 'E': 1, 'T1': 3, 'T2': 0}, 'Odd': {'A1': 5, 'A2': 0, 'E': 1, 'T1': 0, 'T2': 3}}
    print('spin under Td(M):', plain(spin_td), '| physical states by parity:', plain(weights))

    # ---- vibrations: the 15 displacements, then the 9 normal ones (without the 3 translations and the 3 rotations)
    gamma_3n = td(lambda p: chi_3n[p])
    ok &= gamma_3n == {'A1': 1, 'A2': 0, 'E': 1, 'T1': 1, 'T2': 3}
    D5 = [list(d) + [0] for d in zip(*D)]                                   # the directions as a 3 x 5 matrix, carbon at 0
    gens = [[[0, 0, 0], [0, 0, -1], [0, 1, 0]], [[0, 0, 1], [0, 0, 0], [-1, 0, 0]], [[0, -1, 0], [1, 0, 0], [0, 0, 0]]]
    basis = [[1 if r == a else 0 for c in range(5) for r in range(3)] for a in range(3)]
    basis += [[M[r][c] for c in range(5) for r in range(3)] for M in (mat(g, D5) for g in gens)]
    B = T_(basis)
    ok &= rank(B) == 6
    chi_6 = {p: sum(solve(B, [r[0] for r in mat(A15[p], [[x] for x in basis[k]])])[k] for k in range(6)) for p in allp}
    gamma_vib = td(lambda p: chi_3n[p] - chi_6[p])
    ok &= gamma_vib == {'A1': 1, 'A2': 0, 'E': 1, 'T1': 0, 'T2': 2}
    print('displacements under Td(M):', plain(gamma_3n), '| the nine normal ones:', plain(gamma_vib))

    # ---- the J ladder: chi_J(h) = sum_{m=-J}^{J} cos(m theta_h), by the Chebyshev recursion from cos theta_h
    def chi_J(J, c):
        cm = [Fr(1), c]
        while len(cm) <= J: cm.append(2 * c * cm[-1] - cm[-2])
        return cm[0] + 2 * sum(cm[1:J + 1])
    j_table = []
    for J in range(7):
        mult = td(lambda p: chi_J(J, cosang[p]))
        ok &= sum(m * TD[k][(1, 1, 1, 1)] for k, m in mult.items()) == 2 * J + 1
        j_table.append({'J': J, 'species': {k: int(m) for k, m in mult.items() if m},
                        'even': int(sum(m * weights['Even'][k] for k, m in mult.items())), 'odd': int(sum(m * weights['Odd'][k] for k, m in mult.items()))})
    ok &= [r['species'] for r in j_table[:5]] == [{'A1': 1}, {'T1': 1}, {'E': 1, 'T2': 1}, {'A2': 1, 'T1': 1, 'T2': 1}, {'A1': 1, 'E': 1, 'T1': 1, 'T2': 1}]
    print('J ladder:', [(r['J'], r['species'], r['even'], r['odd']) for r in j_table[:5]])

    # ---- T = the even relabelings; 1E(h) = e^{i theta_h}, theta_h the signed turn about the bond of the fixed hydrogen
    def oneE(p):
        if cycles(p) != (3, 1): return W(1)
        f = next(i for i in range(4) if p[i] == i)
        return OMEGA if sum(a * b for a, b in zip(vec(R[p]), D[f])) > 0 else OMEGA * OMEGA
    TT = {'A': lambda p: W(1), '1E': oneE, '2E': lambda p: oneE(p).conj(), 'T': lambda p: W(tr(R[p]))}
    def tt(f): return {k: (sum((W(f(p)) * g(p).conj() for p in even), W(0)) * W(Fr(1, 12))).rational() for k, g in TT.items()}
    ok &= all(oneE(a) * oneE(b) == oneE(comp(a, b)) for a in even for b in even)
    g = next(p for p in even if p[0] == 0 and cycles(p) == (3, 1) and oneE(p) == OMEGA)
    spin_t = tt(lambda p: chi_spin[p])
    ok &= spin_t == {'A': 5, '1E': 1, '2E': 1, 'T': 3}
    reg = tt(lambda p: 12 if p == e else 0)
    ok &= reg == {'A': 1, '1E': 1, '2E': 1, 'T': 3}
    on_t = {k: {n: m for n, m in tt(lambda p, v=v: v[cycles(p)]).items() if m} for k, v in TD.items()}
    ok &= on_t == {'A1': {'A': 1}, 'A2': {'A': 1}, 'E': {'1E': 1, '2E': 1}, 'T1': {'T': 1}, 'T2': {'T': 1}}
    ring = {K: [k for k in ('A', '1E', '2E') if TT[k](g) == [W(1), OMEGA, OMEGA * OMEGA][K % 3]] for K in range(-4, 5)}
    ok &= all(ring[K] == [['A'], ['1E'], ['2E']][K % 3] for K in ring)
    print('1E a homomorphism with 1E(g) = omega | spin under T:', plain(spin_t), '| C[T]:', plain(reg), '| on T:', plain(on_t))

    # ---- isomers, entanglement, monodromy
    iso = []
    for r, fr in TT.items():
        for nu, fn in TT.items():
            if (sum((fr(p) * fn(p) for p in even), W(0)) * W(Fr(1, 12))).rational():
                iso.append((r, nu, int(fr(e).rational()), int(spin_t[nu])))
    swap = {'A': 'A', '1E': '2E', '2E': '1E', 'T': 'T'}
    ok &= iso == [('A', 'A', 1, 5), ('1E', '2E', 1, 1), ('2E', '1E', 1, 1), ('T', 'T', 3, 3)]
    ok &= sorted((swap[r], swap[n], d, m) for r, n, d, m in iso) == sorted(iso)
    kron = lambda A: [[A[i // 3][k // 3] * A[i % 3][k % 3] for k in range(9)] for i in range(9)]
    inv = nullspace([row for p in even for row in [[x - (i == k) for k, x in enumerate(r)] for i, r in enumerate(kron(R[p]))]])
    schmidt = rank([inv[0][3 * i:3 * i + 3] for i in range(3)]) if len(inv) == 1 else None
    entangled = [schmidt * int(spin_t['T']), chi_spin[e]] if schmidt else None
    ok &= schmidt == 3 and entangled == [9, 16]
    mono = {k: 12 // sum(1 for p in even if f(p) == f(e)) for k, f in TT.items()}
    ok &= mono == {'A': 1, '1E': 3, '2E': 3, 'T': 12}
    loop = [next(j for j in range(4) if mat(R[g], [[x] for x in D[i]]) == [[x] for x in D[j]]) for i in range(4)]
    cyclic = sorted(map(sorted, (map(abs, r) for r in R[g]))) == [[0, 0, 1]] * 3 and all(R[g][i][i] == 0 for i in range(3))
    ok &= loop[0] == 0 and sorted(loop) == [0, 1, 2, 3] and loop != [0, 1, 2, 3] and cyclic
    print('isomers:', iso, '| Schmidt rank', schmidt, '| monodromy:', mono, '| the loop returns X0 with', [w + 1 for w in loop])

    # ---- the ball in Rodrigues coordinates r = tan(theta/2) n = vec(R - R^T)/(1 + tr R)
    rod = lambda M: tuple(Fr(x, 1 + tr(M)) for x in vec(M))
    cube = sorted(rod(R[p]) for p in even if cycles(p) == (3, 1))
    ok &= cube == sorted(itertools.product((1, -1), repeat=3))
    # nearer to the identity than to R_h: |<q, 1>| >= |<q, q_h>| with q ~ (1, r), so |1 + r.r_h| <= sqrt(1 + |r_h|^2) for a
    # third-turn (sqrt(1 + 3) = 2) and |r.n| <= 1 for a half-turn about n: r.s <= 1 for every sign vector s, which is
    # |x| + |y| + |z| <= 1, and the rest hold at its six vertices +-e_i, hence on all of it
    verts = [tuple(sgn * (i == k) for k in range(3)) for i in range(3) for sgn in (1, -1)]
    dot = lambda u, v: sum(a * b for a, b in zip(u, v))
    ok &= all(1 + dot(c, c) == 4 for c in cube)
    ok &= all(abs(1 + dot(v, c)) <= 2 for v in verts for c in cube)
    halves = [tuple(Fr(x, max(map(abs, axis(p)))) for x in axis(p)) for p in even if cosang[p] == -1]
    ok &= sorted(map(lambda n: tuple(map(abs, n)), halves)) == [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
    ok &= all(abs(dot(v, n)) <= 1 for v in verts for n in halves)
    quarter = sorted(rod(R[p]) for p in allp if cycles(p) == (4,))
    centers = sorted(tuple(Fr(sum(c[k] for c in cube if c[i] == sg), 4) for k in range(3)) for i in range(3) for sg in (1, -1))
    ok &= quarter == sorted(verts) and centers == sorted(verts)
    def from_rod(r):
        n2 = dot(r, r); K = [[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]]; K2 = mat(K, K)
        return [[(i == k) + Fr(2, 1 + n2) * (K[i][k] + K2[i][k]) for k in range(3)] for i in range(3)]
    Rg = R[g]; glued = True
    for a, b in ((Fr(i, 7), Fr(j, 7)) for i in range(8) for j in range(8 - i)):
        M = from_rod((a, b, 1 - a - b)); Mg = mat(M, T_(Rg))
        glued &= tr(M) == tr(mat(T_(Rg), M)) and tr(Mg) == tr(M) and sum(rod(Mg)) == -1
    twist = rod(mat(from_rod((1, 0, 0)), T_(Rg)))
    ok &= glued and rod(Rg) == (1, 1, 1) and dot(twist, twist) == 1 and twist[0] == 0
    print('ball: the cell is the octahedron |x|+|y|+|z| <= 1, vertices the quarter-turns = the face centers of the cube of',
          'third-turns; the face x+y+z = 1 is glued to x+y+z = -1, the vertex (1,0,0) landing at', tuple(map(int, twist)))

    # ---- packet shares as polynomials in the class overlaps (1, c3, c2)
    def share(k):
        f = TT[k]; out = []
        for cls in ((1, 1, 1, 1), (3, 1), (2, 2)):
            out.append((sum((f(p).conj() for p in even if cycles(p) == cls), W(0)) * f(e) * W(Fr(1, 12))).rational())
        return tuple(out)
    ok &= [share(k) for k in TT] == [(Fr(1, 12), Fr(8, 12), Fr(3, 12)), (Fr(1, 12), Fr(-4, 12), Fr(3, 12)),
                                     (Fr(1, 12), Fr(-4, 12), Fr(3, 12)), (Fr(9, 12), Fr(0), Fr(-9, 12))]
    print('packet shares (coefficients of 1, c3, c2):', plain([share(k) for k in TT]))

    # ---- against the emitted numbers
    S = json.load(open(os.path.join(ROOT, 'figures', 'data', 'summary.json')))
    S0, S1 = S['methane_loop'], S['methane_rigid']
    as_int = lambda d: {k: int(v) for k, v in d.items()}
    ok &= S0['nucleus_i_sits_where_j_was'] == [w + 1 for w in loop] and S0['axis'] == 'C-H1'
    ok &= S1['spin_Td'] == as_int(spin_td) and S1['weights'] == {k: as_int(v) for k, v in weights.items()} and S1['spin_T'] == as_int(spin_t)
    ok &= [tuple(x) for x in S1['isomers']] == iso and S1['entangled_fraction'] == entangled and S1['monodromy_orders'] == mono
    ok &= S1['gamma_3N'] == as_int(gamma_3n) and S1['gamma_vib'] == as_int(gamma_vib) and S1['J_table'] == j_table
    print('summary.json agrees:', ok)
    if not ok:
        sys.exit('verify_methane: mismatch')
    print('verify_methane: all checks pass')

if __name__ == '__main__':
    main()
