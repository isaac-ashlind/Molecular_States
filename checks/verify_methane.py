"""The closing page: rigid methane against Albert, Kubischta, Lemeshko and Liu (arXiv:2403.04572v4), recomputed
from explicit matrices and compared with the numbers make_data emits (summary.json, key methane_rigid).

  geometry   every element of Td(M) = {even relabellings} u {odd relabellings}* carries the tetrahedron X0 into its
             rotational orbit (a proper rotation R_h, residual ~ 0); the twelve unstarred ones form a group of order 12
  spin       (C^2)^4 under the relabellings: 5 A1 + E + 3 T2 (Td labels); under the proper rotations: 5 A + 1E + 2E + 3 T
  pairing    physical states per rotational species and parity (chi_stat chi_pm is A2 for even parity, A1 for odd)
  isomers    Gamma_rot x Gamma_nuc contains A: (A, A) 5, (1E, 2E) 1, (2E, 1E) 1, (T, T) 3, with d = 1, 1, 1, 3
  entangled  the invariant vector of T x T (T = the rotation matrices R_h) has Schmidt rank 3; 9 of the 16 spin states
  monodromy  |T / ker Gamma| = 1, 3, 3, 12
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

def main():
    ok = True
    X0 = G.methane()
    perms = [list(p) for p in itertools.permutations(range(4))]
    # the stabilizer: sigma unstarred when even, starred when odd; R_h from the Kabsch fit on all five nuclei
    # the fit is an exact orthogonal Procrustes by SVD (the suite's Newton polar routine is tuned to small turns and
    # half-turns and does not converge for the third-turns here)
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
    closed = all(any(np.allclose(R[a] @ R[b], R[c]) for c in even) for a in even for b in even)
    ok &= closed and len(even) == 12
    print('geometry: 24 relabellings carry X0 into its orbit by proper rotations; the 12 even ones close:', closed)
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
    ok &= hom
    TT = {'A': lambda q: 1, '1E': oneE, '2E': lambda q: oneE(q).conjugate() if isinstance(oneE(q), complex) else oneE(q),
          'T': lambda q: np.trace(R[q])}
    spin_t = {k: round((sum(chi[q] * complex(f(q)).conjugate() for q in even) / 12).real) for k, f in TT.items()}
    ok &= spin_t == {'A': 5, '1E': 1, '2E': 1, 'T': 3}
    print('1E is a homomorphism of the even relabellings:', hom, '| spin under T:', spin_t)
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
    S = json.load(open(os.path.join(ROOT, 'figures', 'data', 'summary.json')))['methane_rigid']
    ok &= S['spin_Td'] == spin_td and S['weights'] == weights and S['spin_T'] == spin_t
    ok &= [tuple(x) for x in S['isomers']] == iso and S['entangled_fraction'] == [9, 16] and S['monodromy_orders'] == mono
    print('summary.json agrees:', ok)
    if not ok:
        sys.exit('verify_methane: mismatch')
    print('verify_methane: all checks pass')

if __name__ == '__main__':
    main()
