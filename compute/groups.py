"""Permutation-inversion groups for the methylamine and KRb figures.

An element is (perm, star): perm is a tuple of 0-based images, star is 0/1.
The product is composition, (s t)(i) = s(t(i)), which matches the manuscript's
left action s.(t.X) = (st).X with s.X = X P_s and P_{st} = P_t P_s.
"""
import math
from itertools import permutations

def E(n):
    return (tuple(range(n)), 0)

def mul(a, b):
    return (tuple(a[0][b[0][i]] for i in range(len(a[0]))), a[1] ^ b[1])

def inv(a):
    p = [0]*len(a[0])
    for i, j in enumerate(a[0]):
        p[j] = i
    return (tuple(p), a[1])

def cyc(n, *cycles, star=0):
    """Element from 1-based cycles, e.g. cyc(7, (1,2,3)) or cyc(7,(2,3),(4,5),star=1)."""
    p = list(range(n))
    for c in cycles:
        for i, j in zip(c, c[1:] + c[:1]):
            p[i-1] = j-1
    return (tuple(p), star)

def close(gens, n):
    s = {E(n)}; todo = [E(n)]
    while todo:
        a = todo.pop()
        for b in gens:
            c = mul(a, b)
            if c not in s:
                s.add(c); todo.append(c)
    return frozenset(s)

def cosets(G, H):
    """Left cosets gH as frozensets, in a deterministic order."""
    seen = set(); out = []
    for g in sorted(G):
        c = frozenset(mul(g, h) for h in H)
        if c not in seen:
            seen.add(c); out.append(c)
    return out

def interval(H, B, n):
    """All subgroups K with H <= K <= B (brute force by generator extension)."""
    found = {H}; todo = [H]
    while todo:
        k = todo.pop()
        for g in B - k:
            z = close(list(k) + [g], n)
            if z not in found:
                found.add(z); todo.append(z)
    return found

def covers(subs):
    subs = list(subs)
    return [(a, b) for a in subs for b in subs
            if a < b and not any(a < k < b for k in subs)]

def cycle_count(perm):
    seen = set(); c = 0
    for i in range(len(perm)):
        if i not in seen:
            c += 1; j = i
            while j not in seen:
                seen.add(j); j = perm[j]
    return c

def sign(perm):
    return (-1)**(len(perm) - cycle_count(perm))

def cycle_notation(g, n, labels=None):
    """Cycle notation on the first n labels (1-based), with * for inversion."""
    perm, star = g
    seen = set(); parts = []
    for i in range(n):
        if i in seen or perm[i] == i:
            seen.add(i); continue
        c = []; j = i
        while j not in seen:
            seen.add(j); c.append(j); j = perm[j]
        parts.append('(' + ''.join(str(k+1) for k in c) + ')')
    s = ''.join(parts) if parts else 'E'
    return s + ('^*' if star else '')

# ----------------------------------------------------------- methylamine ---

N = 7
b = cyc(N, (2, 3), (4, 5), star=1)
t = cyc(N, (1, 2, 3))
u = cyc(N, (4, 5))
Estar = (tuple(range(N)), 1)
tu = mul(t, u)

H = close([b], N)
G6 = close([b, t], N)
G12 = close([b, t, u], N)
B = close([b, t, u, Estar], N)

def full_S(n_like=5, n=N):
    """S = S_5 x {E, E*} on labels 1..5, fixing 6 and 7."""
    els = set()
    for p in permutations(range(n_like)):
        perm = tuple(p) + tuple(range(n_like, n))
        els.add((perm, 0)); els.add((perm, 1))
    return frozenset(els)

S = full_S()

# Named subgroups of the bond interval [H, B].
NAMED = [
    ('H', [b]),
    ('\\langle H,u\\rangle', [b, u]),
    ('\\langle H,E^*\\rangle', [b, Estar]),
    ('\\langle H,u^*\\rangle', [b, mul(u, Estar)]),
    ('G_6', [b, t]),
    ('\\langle H,u,E^*\\rangle', [b, u, Estar]),
    ('G_{12}', [b, t, u]),
    ('\\langle G_6,E^*\\rangle', [b, t, Estar]),
    ('\\langle G_6,u^*\\rangle', [b, t, mul(u, Estar)]),
    ('B', [b, t, u, Estar]),
]

def interval_data():
    subs = interval(H, B, N)
    names = {}
    for name, gens in NAMED:
        K = close(gens, N)
        assert K in subs, name
        names[K] = name
    assert len(names) == len(subs) == 10
    cov = covers(subs)
    return subs, names, cov

# G6 character table on classes {E}, {t, t^2}, {b, tb, t^2b}
def g6_class(g):
    if g == E(N):
        return 'E'
    if g[1] == 0:
        return 't'
    return 'b'

G6_CHARS = {'A1': {'E': 1, 't': 1, 'b': 1},
            'A2': {'E': 1, 't': 1, 'b': -1},
            'E':  {'E': 2, 't': -1, 'b': 0}}

def spin_weights(G=G6, n_protons=5, spin=0.5):
    """Multiplicity of each partner spin species, by total parity.

    weight[parity][Gamma] = (1/|G|) sum_g chi_Gamma(g) chi_spin(g) chi_stat(g) chi_parity(g)
    with chi_spin(g) = (2I+1)^{cycles of sigma_g on the protons},
    chi_stat(g) = (sgn sigma_g)^{2I}, chi_parity(g) = parity^{star(g)}.
    It equals the number of copies of Gamma (x) chi in the spin space, i.e. the
    number of physical states per spatial level of species Gamma.
    """
    out = {}
    for parity in (1, -1):
        vals = {}
        for name, chi in G6_CHARS.items():
            tot = 0
            for g in G:
                perm = g[0][:n_protons]
                chi_spin = (int(round(2*spin))+1)**cycle_count(perm)
                chi_stat = sign(perm)**int(round(2*spin))
                chi_par = parity**g[1]
                tot += chi[g6_class(g)]*chi_spin*chi_stat*chi_par
            assert tot % len(G) == 0
            vals[name] = tot // len(G)
        out[parity] = vals
    return out

def coset_action_table():
    """Left action of t and b on the ordered cosets (H, tH, t^2H)."""
    cs = [H, frozenset(mul(t, h) for h in H), frozenset(mul(mul(t, t), h) for h in H)]
    names = ['H', 'tH', 't^2H']
    def act(g):
        return [cs.index(frozenset(mul(g, x) for x in c)) for c in cs]
    return names, act(t), act(b)

# ------------------------------------------------------------------- KRb ---

def krb_groups():
    n = 4
    p12 = cyc(n, (1, 2)); p34 = cyc(n, (3, 4)); Es = (tuple(range(n)), 1)
    Sk = close([p12, p34, Es], n)
    Gin = close([mul(p12, p34), Es], n)
    GK = close([p12, Es], n)
    GRb = close([p34, Es], n)
    return Sk, Gin, GK, GRb, p12, p34, Es
