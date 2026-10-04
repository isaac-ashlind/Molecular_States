"""Centered nuclear configurations, group actions on them, and projections.

Everything here is plain Python (no numpy).  A configuration X is stored as a
list of N columns, each a 3-vector (x, y, z) in angstrom.  Masses are in
unified atomic mass units.  The conventions follow the manuscript:

    s . X = X P_s,     P_{sigma*} = -P_sigma,     (X P_sigma)_i = x_{sigma^{-1}(i)},
    r . X = R X,       P_{st} = P_t P_s.

Labels are 1-based in the public API, matching the manuscript's labels.
"""
import math

# Isotope masses (u).  Sources: AME2020 rounded to 5 decimals.
MASS = {
    'H': 1.00783,    # 1H
    'C': 12.00000,   # 12C
    'N': 14.00307,   # 14N
    'O': 15.99491,   # 16O
    'K': 39.96400,   # 40K
    'Rb': 86.90918,  # 87Rb
}

# Nuclear spin numbers used in the figures.
SPIN = {'H': 0.5, 'C': 0.0, 'N': 1.0, 'O': 0.0, 'K': 4.0, 'Rb': 1.5}

# ---------------------------------------------------------------- vectors ---

def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def scale(a, s): return tuple(x * s for x in a)
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def norm(a): return math.sqrt(dot(a, a))
def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def unit(a):
    n = norm(a)
    return tuple(x / n for x in a)

def matvec(M, v):
    return tuple(dot(row, v) for row in M)

def matmul(A, B):
    n = len(A); m = len(B[0]); k = len(B)
    return [[sum(A[i][l]*B[l][j] for l in range(k)) for j in range(m)] for i in range(n)]

def transpose(A):
    return [list(col) for col in zip(*A)]

def identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def rot_x(a):
    c, s = math.cos(a), math.sin(a)
    return [[1, 0, 0], [0, c, -s], [0, s, c]]

def rot_y(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, 0, s], [0, 1, 0], [-s, 0, c]]

def rot_z(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]]

def rot_axis(axis, a):
    """Rodrigues rotation about a unit axis by angle a."""
    x, y, z = unit(axis)
    c, s = math.cos(a), math.sin(a)
    C = 1 - c
    return [[c + x*x*C, x*y*C - z*s, x*z*C + y*s],
            [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
            [z*x*C - y*s, z*y*C + x*s, c + z*z*C]]

# ---------------------------------------------------------- configurations ---

def center(X, masses):
    """Translate so that X m = 0 (mass-weighted centre at the origin)."""
    M = sum(masses)
    c = (sum(m*x[0] for m, x in zip(masses, X)) / M,
         sum(m*x[1] for m, x in zip(masses, X)) / M,
         sum(m*x[2] for m, x in zip(masses, X)) / M)
    return [sub(x, c) for x in X]

def mass_moment(X, masses):
    """The 3-vector X m; zero for a centred configuration."""
    return tuple(sum(m*x[k] for m, x in zip(masses, X)) for k in range(3))

def rotate(R, X):
    return [matvec(R, x) for x in X]

def invert_config(X):
    return [scale(x, -1.0) for x in X]

def perm_matrix(sigma, n):
    """Column-permutation matrix P_sigma with (X P)_i = x_{sigma^{-1}(i)}.

    sigma is a dict or list mapping 0-based i -> sigma(i).  P[j][i] = 1 iff
    j = sigma^{-1}(i), i.e. sigma(j) = i.
    """
    P = [[0]*n for _ in range(n)]
    for j in range(n):
        P[j][sigma[j]] = 1
    return P

def apply_perm_inversion(X, sigma, star):
    """s . X = X P_s with P_{sigma*} = -P_sigma.  sigma is 0-based list."""
    n = len(X)
    inv = [0]*n
    for j in range(n):
        inv[sigma[j]] = j
    Y = [X[inv[i]] for i in range(n)]
    if star:
        Y = invert_config(Y)
    return Y

def config_distance(X, Y, masses=None):
    """Euclidean (or mass-weighted) distance between two configurations."""
    if masses is None:
        masses = [1.0]*len(X)
    return math.sqrt(sum(m*dot(sub(x, y), sub(x, y)) for m, x, y in zip(masses, X, Y)))

def mass_inner(A, B, masses):
    """<A,B>_m = sum_i m_i a_i . b_i for displacement patterns A, B."""
    return sum(m*dot(a, b) for m, a, b in zip(masses, A, B))

def kabsch_rotation(X, Y, masses):
    """Return (R, residual): best proper rotation with Y ~ R X in mass metric.

    Implemented by a small Jacobi SVD-free method: we solve the orthogonal
    Procrustes problem via the polar decomposition of the 3x3 cross-covariance
    computed with Newton iteration for the symmetric square root.  Adequate for
    the well-conditioned, exactly-symmetric cases used here (residual ~1e-12).
    """
    H = [[sum(m*x[i]*y[j] for m, x, y in zip(masses, X, Y)) for j in range(3)] for i in range(3)]
    # polar decomposition H = U S, U orthogonal, by Higham iteration
    U = [row[:] for row in H]
    for _ in range(60):
        Ui = inverse3(U)
        UiT = transpose(Ui)
        U = [[0.5*(U[i][j] + UiT[i][j]) for j in range(3)] for i in range(3)]
    # U is orthogonal with U = H (H^T H)^{-1/2}; Y ~ R X needs R = U
    R = U
    if det3(R) < 0:
        # improper: not a rotation; reflect the last column so the caller sees a large residual
        R = [[R[i][j] * (-1 if j == 2 else 1) for j in range(3)] for i in range(3)]
    res = config_distance(rotate(R, X), Y, masses)
    return R, res

def det3(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
            - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
            + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

def inverse3(M):
    d = det3(M)
    a, b, c = M[0]; dd, e, f = M[1]; g, h, i = M[2]
    return [[(e*i-f*h)/d, (c*h-b*i)/d, (b*f-c*e)/d],
            [(f*g-dd*i)/d, (a*i-c*g)/d, (c*dd-a*f)/d],
            [(dd*h-e*g)/d, (b*g-a*h)/d, (a*e-b*dd)/d]]

# ------------------------------------------------------------------- water ---

WATER_ELEMENTS = ['H', 'H', 'O']            # labels 1, 2, 3
WATER_MASSES = [MASS[e] for e in WATER_ELEMENTS]
WATER_ELL = 0.9575                           # mean O-H bond length (angstrom)

def water(delta_ratio, theta_deg, ell=WATER_ELL):
    """Centred planar water with r1 = ell(1+delta_ratio), r2 = ell(1-delta_ratio).

    Oxygen is nucleus 3.  Before centring the oxygen sits at the origin and the
    bisector of the bond angle points along -y, so the hydrogens are above.
    The returned configuration is centred on the mass centre, which is NOT the
    oxygen position: that is the point of drawing it this way.
    """
    th = math.radians(theta_deg)
    r1 = ell*(1+delta_ratio); r2 = ell*(1-delta_ratio)
    H1 = (-r1*math.sin(th/2), r1*math.cos(th/2), 0.0)
    H2 = ( r2*math.sin(th/2), r2*math.cos(th/2), 0.0)
    O = (0.0, 0.0, 0.0)
    return center([H1, H2, O], WATER_MASSES)

# ------------------------------------------------------------ methylamine ---

MLA_ELEMENTS = ['H', 'H', 'H', 'H', 'H', 'C', 'N']   # labels 1..7
MLA_MASSES = [MASS[e] for e in MLA_ELEMENTS]
MLA_BONDS = [(6, 7), (6, 1), (6, 2), (6, 3), (7, 4), (7, 5)]   # 1-based

# Model parameters for the reference family (angstrom, degrees).  These are
# representative structural values for CH3NH2, not a fitted equilibrium
# geometry.  The family X0(tau, eta) follows the manuscript's footnote:
# methyl hydrogens at angles tau - 2 pi (k-1)/3 at fixed radius and height
# about the C-N axis, amino hydrogens at (eta, +-a, z_A).
MLA = {
    'r_CN': 1.471, 'r_CH': 1.093, 'r_NH': 1.010,
    'angle_HCN': 110.3, 'angle_HNC': 110.0, 'angle_HNH': 107.0,
}

def _mla_frame_params(p=MLA):
    rho_M = p['r_CH']*math.sin(math.radians(p['angle_HCN']))
    z_M = p['r_CH']*math.cos(math.radians(p['angle_HCN']))      # below carbon (negative)
    z_C = 0.0
    z_N = p['r_CN']
    proj = p['r_NH']*math.sin(math.radians(p['angle_HNC']))
    z_A = z_N + p['r_NH']*(-math.cos(math.radians(p['angle_HNC'])))   # above nitrogen
    dz = z_A - z_N
    # projected half-angle phi between the two N-H bonds: solve the HNH angle
    cos_hnh = math.cos(math.radians(p['angle_HNH']))
    cos2phi = (cos_hnh*p['r_NH']**2 - dz*dz)/(proj*proj)
    phi = 0.5*math.acos(cos2phi)
    eta0 = proj*math.cos(phi)
    a = proj*math.sin(phi)
    return dict(rho_M=rho_M, z_M=z_M, z_C=z_C, z_N=z_N, z_A=z_A, eta0=eta0, a=a, phi=phi)

MLA_FRAME = _mla_frame_params()

def methylamine(tau, eta, frame=MLA_FRAME):
    """Centred reference family X0(tau, eta); tau in radians, eta in angstrom.

    Column order is the manuscript's: H1 H2 H3 (methyl), H4 H5 (amino), C6, N7.
    At tau = 0 the nucleus H1 lies in the +x half of the xz-plane, which is the
    reference mirror plane through H1, C6, N7.
    """
    f = frame
    cols = []
    for k in range(3):
        phi = tau - 2*math.pi*k/3
        cols.append((f['rho_M']*math.cos(phi), f['rho_M']*math.sin(phi), f['z_M']))
    cols.append((eta, f['a'], f['z_A']))
    cols.append((eta, -f['a'], f['z_A']))
    cols.append((0.0, 0.0, f['z_C']))
    cols.append((0.0, 0.0, f['z_N']))
    return center(cols, MLA_MASSES)

def methylamine_uncentred(tau, eta, frame=MLA_FRAME):
    f = frame
    cols = []
    for k in range(3):
        phi = tau - 2*math.pi*k/3
        cols.append((f['rho_M']*math.cos(phi), f['rho_M']*math.sin(phi), f['z_M']))
    cols.append((eta, f['a'], f['z_A']))
    cols.append((eta, -f['a'], f['z_A']))
    cols.append((0.0, 0.0, f['z_C']))
    cols.append((0.0, 0.0, f['z_N']))
    return cols

def d_methylamine(tau, eta, h=1e-5):
    """Tangent vectors (dX0/dtau, dX0/deta) by central differences."""
    Xp = methylamine(tau+h, eta); Xm = methylamine(tau-h, eta)
    dt = [scale(sub(a, b), 1/(2*h)) for a, b in zip(Xp, Xm)]
    Xp = methylamine(tau, eta+h); Xm = methylamine(tau, eta-h)
    de = [scale(sub(a, b), 1/(2*h)) for a, b in zip(Xp, Xm)]
    return dt, de

def rotation_generators(X):
    """Infinitesimal rotations L_k X for k = x, y, z (displacement patterns)."""
    axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    return [[cross(ax, x) for x in X] for ax in axes]

def gram_schmidt_mass(vectors, masses, tol=1e-9):
    """Mass-orthonormalise a list of displacement patterns (drops dependent ones)."""
    out = []
    for v in vectors:
        w = [tuple(x) for x in v]
        for u in out:
            c = mass_inner(w, u, masses)
            w = [sub(a, scale(b, c)) for a, b in zip(w, u)]
        n = math.sqrt(mass_inner(w, w, masses))
        if n > tol:
            out.append([scale(a, 1/n) for a in w])
    return out

def normal_vector(pattern, tau, eta, masses=MLA_MASSES):
    """Project a displacement pattern onto the mass-normal space at X0(tau,eta).

    The normal space is the mass-orthogonal complement of the shape tangents
    and the three rotation generators; the pattern is first made translation
    free (sum m_i v_i = 0).  Returns the unit normal vector.
    """
    X = methylamine(tau, eta)
    dt, de = d_methylamine(tau, eta)
    tangent = gram_schmidt_mass([dt, de] + rotation_generators(X), masses)
    M = sum(masses)
    c = scale(tuple(sum(m*v[k] for m, v in zip(masses, pattern)) for k in range(3)), 1/M)
    w = [sub(v, c) for v in pattern]
    for u in tangent:
        cc = mass_inner(w, u, masses)
        w = [sub(a, scale(b, cc)) for a, b in zip(w, u)]
    n = math.sqrt(mass_inner(w, w, masses))
    return [scale(a, 1/n) for a in w], tangent

# ------------------------------------------------------------------- KRb ---

KRB_ELEMENTS = ['K', 'K', 'Rb', 'Rb']   # labels 1, 2, 3, 4
KRB_MASSES = [MASS[e] for e in KRB_ELEMENTS]

def krb_schematic():
    """A schematic planar K2Rb2 arrangement (angstrom-like units), centred.

    This is NOT a computed complex geometry; it is declared schematic in the
    figure.  Potassium on the upper row, rubidium on the lower row, so that the
    three groupings (incoming pairs K1Rb3 | K2Rb4, one complex, outgoing pairs
    K1K2 | Rb3Rb4) are drawn by envelopes that never cross.
    """
    cols = [(-1.35, 0.8, 0.0),   # K1
            ( 1.35, 0.8, 0.0),   # K2
            (-1.35, -0.8, 0.0),  # Rb3
            ( 1.35, -0.8, 0.0)]  # Rb4
    return center(cols, KRB_MASSES)

# ------------------------------------------------------------- projection ---

def camera(azimuth_deg, elevation_deg):
    """Orthographic camera.  Returns (right, up, towards-viewer) unit vectors.

    The view direction is rotated from +y (looking at the xz-plane face-on,
    with x up the page and z to the right) by azimuth about z and elevation
    towards +x ... concretely: page_x = molecular z, page_y = x cos b + y sin b
    for elevation 0 and azimuth b; elevation tilts the C-N axis.
    """
    b = math.radians(azimuth_deg); e = math.radians(elevation_deg)
    # basis before elevation: right = z, up = (cos b, sin b, 0), out = up x right
    right = (0.0, 0.0, 1.0)
    up = (math.cos(b), math.sin(b), 0.0)
    out = cross(up, right)   # towards the viewer
    # elevation: rotate right and out about the up axis
    R = rot_axis(up, e)
    right = matvec(R, right); out = matvec(R, out)
    return right, up, out

def project(X, cam):
    right, up, out = cam
    return [(dot(x, right), dot(x, up), dot(x, out)) for x in X]   # (px, py, depth)

def draw_order(P):
    """Indices sorted back to front (increasing depth = towards viewer)."""
    return sorted(range(len(P)), key=lambda i: P[i][2])
