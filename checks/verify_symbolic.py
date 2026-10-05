"""Symbolic (exact) checks of the identities behind the figures.  Needs sympy.

Run from the repository root:  python3 checks/verify_symbolic.py
The actions of t, u, b, tu on the family X0(tau, iota) hold identically; the overlap c = exp(-d^2/8 Delta^2) and the
density variance Delta^2; the shares (1 + 2c)/3 = 1/3 + 2c/3 and 2(1 - c)/3 from the projectors on three packets of
overlap c, with P_A1 eta_0 the average of the three; two neighboring densities cross at c times the peak;
the E-pair matrices of t and b; for a packet with arbitrary overlaps, the shares of its parts even and odd under b
add, the odd part carries only A2 and E with weight (1 - <eta, U_b eta>)/2, and with overlaps zero outside G6 each
G6 share splits evenly over the two G12 species above it; the seam relation kappa = m + rho K for every integer m.
"""
import sympy as sp

tau, iota, rhoM, zM, bA, zA, zC, zN, mH, mC, mN = sp.symbols('tau iota rho_M z_M b_A z_A z_C z_N m_H m_C m_N', real=True)
masses = [mH]*5 + [mC, mN]

def family(tau, iota):
    """The reference family before centering (centering commutes with the actions)."""
    cols = []
    for k in range(3):
        phi = tau - 2*sp.pi*k/3
        cols.append(sp.Matrix([rhoM*sp.cos(phi), rhoM*sp.sin(phi), zM]))
    cols += [sp.Matrix([iota, bA, zA]), sp.Matrix([iota, -bA, zA]), sp.Matrix([0, 0, zC]), sp.Matrix([0, 0, zN])]
    return cols

def center(cols):
    M = sum(masses)
    c = sum((m*x for m, x in zip(masses, cols)), sp.zeros(3, 1))/M
    return [x - c for x in cols]

def act(cols, perm, star):
    """s.X = X P_s: column i receives x_{sigma^{-1}(i)}; starred: negate."""
    n = len(cols); inv = [0]*n
    for j in range(n):
        inv[perm[j]] = j
    out = [cols[inv[i]] for i in range(n)]
    return [-x for x in out] if star else out

def R_y(th): return sp.Matrix([[sp.cos(th), 0, sp.sin(th)], [0, 1, 0], [-sp.sin(th), 0, sp.cos(th)]])
def R_z(th): return sp.Matrix([[sp.cos(th), -sp.sin(th), 0], [sp.sin(th), sp.cos(th), 0], [0, 0, 1]])

t = ([1, 2, 0, 3, 4, 5, 6], 0)          # (123): 1->2->3->1 on 0-based labels
u = ([0, 1, 2, 4, 3, 5, 6], 0)          # (45)
b = ([0, 2, 1, 4, 3, 5, 6], 1)          # (23)(45)*
def compose(s, tt):                       # (s t)(i) = s(t(i))
    return ([s[0][tt[0][i]] for i in range(7)], s[1] ^ tt[1])
tu = compose(t, u)

X = center(family(tau, iota))
cases = {
    't':  (t,  (tau + 2*sp.pi/3, iota), sp.eye(3)),
    'u':  (u,  (tau + sp.pi, -iota), R_z(sp.pi)),
    'b':  (b,  (-tau, iota), R_y(sp.pi)),
    'tu': (tu, (tau - sp.pi/3, -iota), R_z(sp.pi)),
}
for name, (g, (tp, ip), A) in cases.items():
    lhs = act(X, g[0], g[1])
    rhs = [A*x for x in center(family(tp, ip))]
    for l, r in zip(lhs, rhs):
        d = sp.simplify(sp.expand_trig(l - r))
        assert d == sp.zeros(3, 1), (name, d)
    print(f'{name}: g.X0(tau, iota) = A_g X0(tau\', iota\') holds identically')
# Gaussian overlap: amplitudes with density variance Delta^2 in the plane
x, y, d, sg = sp.symbols('x y d Delta', real=True, positive=True)
g0 = sp.exp(-(x**2 + y**2)/(4*sg**2)); g1 = sp.exp(-((x - d)**2 + y**2)/(4*sg**2))
num = sp.integrate(sp.integrate(g0*g1, (x, -sp.oo, sp.oo)), (y, -sp.oo, sp.oo))
den = sp.integrate(sp.integrate(g0*g0, (x, -sp.oo, sp.oo)), (y, -sp.oo, sp.oo))
c = sp.simplify(num/den)
assert sp.simplify(c - sp.exp(-d**2/(8*sg**2))) == 0
var = sp.integrate(sp.integrate(x**2*g0**2, (x, -sp.oo, sp.oo)), (y, -sp.oo, sp.oo))/den
assert sp.simplify(var - sg**2) == 0
print('Gaussian: density variance Delta^2 and c = exp(-d^2/(8 Delta^2)) hold exactly')
cc = sp.symbols('c', real=True)
# the shares from the projectors on the span of three packets with pairwise overlap c: the Gram matrix of the packets
# and the site action of t; w_Gamma = <eta_0, P_Gamma eta_0> with P_Gamma = (d_Gamma/3) sum_k chi_Gamma(t^k) t^k
Gm = sp.Matrix(3, 3, lambda i, j: 1 if i == j else cc)
T3 = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
e0 = sp.Matrix([1, 0, 0]); w3 = sp.Rational(-1, 2) + sp.sqrt(3)*sp.I/2   # omega = e^{2 pi i/3}
def share(chars, d):
    P = sum((d*sp.conjugate(chars[k])/3*T3**k for k in range(3)), sp.zeros(3, 3))
    return sp.simplify(sp.expand((e0.T*Gm*P*e0)[0]))
wA1 = share([1, 1, 1], 1); wE = share([1, w3, w3**2], 1) + share([1, w3**2, w3], 1)
assert sp.simplify(wA1 - (1 + 2*cc)/3) == 0 and sp.simplify(wE - 2*(1 - cc)/3) == 0 and sp.simplify(wA1 + wE - 1) == 0
# P_A1 averages eta_0 over the three versions, so with no overlap the A1 share is 1/3 and overlap adds 2c/3
PA1 = sum((T3**k for k in range(3)), sp.zeros(3, 3))/3
assert PA1**2 == PA1 and PA1*e0 == sp.Matrix([1, 1, 1])/3
assert sp.simplify((e0.T*Gm*PA1*e0)[0] - wA1) == 0 and (e0.T*Gm.subs(cc, 0)*PA1*e0)[0] == sp.Rational(1, 3)
# the drawn profiles are densities exp(-x^2/(2 Delta^2)); two neighbors d apart cross at the midpoint, at c times the peak
assert sp.simplify(sp.exp(-(d/2)**2/(2*sg**2)) - c) == 0
print('shares: w_A1 = (1+2c)/3 = 1/3 + 2c/3, w_E = 2(1-c)/3, sum rule, P_A1 eta_0 the average of the three, densities cross at c: exact')

# E pair matrices in the orthonormal basis v2, v3
v1 = sp.Matrix([1, 1, 1])/sp.sqrt(3); v2 = sp.Matrix([2, -1, -1])/sp.sqrt(6); v3 = sp.Matrix([0, 1, -1])/sp.sqrt(2)
T = sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]]); Bm = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
V = sp.Matrix.hstack(v2, v3)
assert (V.T*V - sp.eye(2)).applyfunc(sp.simplify) == sp.zeros(2, 2) and sp.simplify(v1.dot(v2)) == 0 and sp.simplify(v1.dot(v3)) == 0
TE = (V.T*T*V).applyfunc(sp.simplify); BE = (V.T*Bm*V).applyfunc(sp.simplify)
assert TE == sp.Matrix([[-sp.Rational(1, 2), -sp.sqrt(3)/2], [sp.sqrt(3)/2, -sp.Rational(1, 2)]])
assert BE == sp.diag(1, -1)
assert T*v1 == v1 and Bm*v1 == v1
print('E pair: t -> rotation by 2pi/3, b -> diag(1,-1); A1 vector fixed: exact')

# a packet with arbitrary real overlaps f(g) = <eta, U_g eta>, f(E) = 1, f(g^-1) = f(g), U_g U_h = U_gh
def close(gens):
    e = (list(range(7)), 0); out = [e]; todo = [e]
    while todo:
        x = todo.pop()
        for g in gens:
            y = compose(x, g)
            if y not in out:
                out.append(y); todo.append(y)
    return out
def inv(g):
    p = [0]*7
    for i, j in enumerate(g[0]):
        p[j] = i
    return (p, g[1])
def order(g):
    k, x = 1, g
    while x != (list(range(7)), 0):
        x, k = compose(x, g), k + 1
    return k
G6 = close([t, b]); G12 = close([t, b, u])
assert len(G6) == 6 and len(G12) == 12 and all(compose(u, g) == compose(g, u) for g in G12)   # G12 = G6 x <u>
key = lambda g: (tuple(g[0]), g[1])
f = {}
for g in G6:
    f.setdefault(min(key(g), key(inv(g))), sp.Symbol('f%d' % len(f)))
F = lambda g: f[min(key(g), key(inv(g)))] if key(g) in {key(x) for x in G6} else 0   # zero outside G6
f[key((list(range(7)), 0))] = 1
CH = {'A1': {1: 1, 3: 1, 2: 1}, 'A2': {1: 1, 3: 1, 2: -1}, 'E': {1: 2, 3: -1, 2: 0}}   # by element order: E, t class, b class
def shares(ov, group=G6, char=lambda n, g: CH[n][order(g)], names=('A1', 'A2', 'E')):
    """||P_Gamma eta||^2 = (d/|G|) sum_g chi(g) <eta, U_g eta> (real characters)."""
    return {n: sp.expand(sp.Rational(char(n, group[0]), len(group))*sum(char(n, g)*ov(g) for g in group)) for n in names}
def part(sg):          # <eta_s, U_g eta_s> for eta_s = (eta + sg U_b eta)/2
    return lambda g: sp.Rational(1, 4)*(F(g) + sg*F(compose(g, b)) + sg*F(compose(b, g)) + F(compose(compose(b, g), b)))
w, wp, wm = shares(F), shares(part(1)), shares(part(-1))
assert all(sp.expand(w[n] - wp[n] - wm[n]) == 0 for n in w) and wp['A2'] == 0 and wm['A1'] == 0
assert sp.expand(sum(wm.values()) - (1 - F(b))/2) == 0
cp = sp.simplify(part(1)(t)/part(1)((list(range(7)), 0)))                 # the even part's own overlap
assert sp.simplify(wp['A1']/sum(wp.values()) - (1 + 2*cp)/3) == 0
print('packet not fixed by H: shares of eta+ and eta- add, eta+ has no A2 and follows (1+2c)/3, eta- has no A1, ||eta-||^2 = (1 - <eta, U_b eta>)/2: exact')
# leakage: with overlaps zero outside G6, the shares under G12 = G6 x <u> follow by branching, half to each sign of u
in6 = {key(x) for x in G6}
ch12 = lambda n, g: CH[n[:-1]][order(g if key(g) in in6 else compose(g, u))] * (1 if key(g) in in6 or n[-1] == '+' else -1)
w12 = shares(F, G12, ch12, [n + sg for n in CH for sg in '+-'])
assert all(sp.expand(w12[n + sg] - w[n]/2) == 0 for n in CH for sg in '+-')
print('overlaps zero outside G6: each G6 share splits evenly between the two G12 species above it: exact')

# seam and Wigner phase: D(r R_z(omega)) = e^{-iK omega} D(r) with (r,tau+2pi,q) ~ (r R_z(-2 pi rho),tau,q)
K, rho, kappa = sp.symbols('K rho kappa'); m = sp.symbols('m', integer=True)
lhs = sp.exp(2*sp.pi*sp.I*kappa)                       # from e^{i kappa (tau + 2 pi)}
rhs = sp.exp(-sp.I*K*(-2*sp.pi*rho))                    # from D(r R_z(-2 pi rho))
assert sp.simplify(sp.powsimp(sp.expand((lhs/rhs).subs(kappa, m + rho*K)))) == 1   # for every integer m
print('seam: e^{2 pi i kappa} = e^{2 pi i rho K} with kappa = m + rho K, m integer: exact')
print('verify_symbolic.py: all checks passed')
