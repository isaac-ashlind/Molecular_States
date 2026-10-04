"""Symbolic (exact) verification of the identities behind the figures.  Needs sympy.

Run from the repository root:  python3 checks/verify_symbolic.py
"""
import sympy as sp

tau, eta, rhoM, zM, a, zA, zC, zN, mH, mC, mN = sp.symbols('tau eta rho_M z_M a z_A z_C z_N m_H m_C m_N', real=True)
masses = [mH]*5 + [mC, mN]

def family(tau, eta):
    """Uncentred reference family; centring commutes with the checks below (see docs)."""
    cols = []
    for k in range(3):
        phi = tau - 2*sp.pi*k/3
        cols.append(sp.Matrix([rhoM*sp.cos(phi), rhoM*sp.sin(phi), zM]))
    cols += [sp.Matrix([eta, a, zA]), sp.Matrix([eta, -a, zA]), sp.Matrix([0, 0, zC]), sp.Matrix([0, 0, zN])]
    return cols

def centre(cols):
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

X = centre(family(tau, eta))
cases = {
    't':  (t,  (tau + 2*sp.pi/3, eta), sp.eye(3)),
    'u':  (u,  (tau + sp.pi, -eta), R_z(sp.pi)),
    'b':  (b,  (-tau, eta), R_y(sp.pi)),
    'tu': (tu, (tau - sp.pi/3, -eta), R_z(sp.pi)),
}
for name, (g, (tp, ep), A) in cases.items():
    lhs = act(X, g[0], g[1])
    rhs = [A*x for x in centre(family(tp, ep))]
    for l, r in zip(lhs, rhs):
        d = sp.simplify(sp.expand_trig(l - r))
        assert d == sp.zeros(3, 1), (name, d)
    print(f'{name}: g.X0(tau,eta) = A_g X0(tau\',eta\') holds identically  (A_g = {"I" if A == sp.eye(3) else name and ("R_z(pi)" if A == R_z(sp.pi) else "R_y(pi)")})')
# centring: X m = 0 identically
Mtot = sum(masses)
assert sp.simplify(sum((m*x for m, x in zip(masses, X)), sp.zeros(3, 1))) == sp.zeros(3, 1)
print('centred family: X m = 0 identically')

# Gaussian overlap: amplitudes with density variance Delta^2 in the plane (width Delta, approved notation)
x, y, d, sg = sp.symbols('x y d Delta', real=True, positive=True)
g0 = sp.exp(-(x**2 + y**2)/(4*sg**2)); g1 = sp.exp(-((x - d)**2 + y**2)/(4*sg**2))
num = sp.integrate(sp.integrate(g0*g1, (x, -sp.oo, sp.oo)), (y, -sp.oo, sp.oo))
den = sp.integrate(sp.integrate(g0*g0, (x, -sp.oo, sp.oo)), (y, -sp.oo, sp.oo))
c = sp.simplify(num/den)
assert sp.simplify(c - sp.exp(-d**2/(8*sg**2))) == 0
var = sp.integrate(sp.integrate(x**2*g0**2, (x, -sp.oo, sp.oo)), (y, -sp.oo, sp.oo))/den
assert sp.simplify(var - sg**2) == 0
print('Gaussian: density variance Delta^2 and c = exp(-d^2/(8 Delta^2)) hold exactly')
cc = sp.symbols('c')
wA1 = (1 + 2*cc)/3; wE = 2*(1 - cc)/3
assert sp.simplify(wA1 + wE - 1) == 0
# shares from the projector with all pairwise overlaps c: ||(g0+g1+g2)/3||^2 = (3 + 6c)/9
assert sp.simplify(sp.Rational(1, 9)*(3 + 6*cc) - wA1) == 0
print('shares: w_A1 = (1+2c)/3, w_E = 2(1-c)/3, sum rule exact')

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

# seam and Wigner phase: D(r R_z(omega)) = e^{-iK omega} D(r) with (r,tau+2pi,q) ~ (r R_z(-2 pi rho),tau,q)
K, m, rho, kappa = sp.symbols('K m rho kappa')
lhs = sp.exp(2*sp.pi*sp.I*kappa)                       # from e^{i kappa (tau + 2 pi)}
rhs = sp.exp(-sp.I*K*(-2*sp.pi*rho))                    # from D(r R_z(-2 pi rho))
assert sp.simplify((lhs - rhs).subs(kappa, m + rho*K).subs(m, 3)) == 0
assert sp.simplify((lhs - rhs).subs(kappa, m + rho*K).subs(m, -2)) == 0
print('seam: e^{2 pi i kappa} = e^{2 pi i rho K} with kappa = m + rho K, m integer: exact')
print('verify_symbolic.py: all checks passed')
