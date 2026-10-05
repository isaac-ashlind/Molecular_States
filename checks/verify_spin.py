"""Spin weights of methylamine (Section 8, figure 8), from explicit permutation matrices on the spin basis.

  the five protons (C^2)^{x5} under G6 = <t = (123), (23)(45)> decompose as 12 A1 + 4 A2 + 8 E, and so does the product
  of the methyl multiplets (C^2)^{x3} under S3 = 4 A1 + 2 E with the amino ones (C^2)^{x2} under (45) = 3 (+) + 1 (-)
  the full spin space (C^2)^{x5} x C^1 x C^3 (carbon-12 and nitrogen-14, not permuted) is 36 A1 + 12 A2 + 24 E
  the physical weight of a spatial level of species Gamma is the multiplicity of the spin species Gamma' whose product
  with Gamma contains chi_stat chi_pm (chi_stat is trivial on G6): 36, 12, 24 for even parity, 12, 36, 24 for odd
Requires numpy. Exit status nonzero on any mismatch.
"""
import sys
import numpy as np

def perm_matrix(perm, n):
    dim = 2 ** n
    M = np.zeros((dim, dim))
    for idx in range(dim):
        bits = [(idx >> k) & 1 for k in range(n)]
        new = [0] * n
        for i in range(n):
            new[perm[i]] = bits[i]
        M[sum(b << k for k, b in enumerate(new)), idx] = 1
    return M

CLS = [1, 2, 3]  # class sizes of S3 (and of G6): identity, the 3-cycles, the transpositions (b and its conjugates)
TAB = {'A1': [1, 1, 1], 'A2': [1, 1, -1], 'E': [2, -1, 0]}

def inner(a, b):
    return round(sum(c * x * y for c, x, y in zip(CLS, a, b)) / 6)

def decompose(chars):
    return {k: inner(chars, v) for k, v in TAB.items()}

def main():
    ok = True
    e3, c3, s3 = {0: 0, 1: 1, 2: 2}, {0: 1, 1: 2, 2: 0}, {0: 0, 1: 2, 2: 1}
    methyl = decompose([np.trace(perm_matrix(p, 3)) for p in (e3, c3, s3)])
    ok &= methyl == {'A1': 4, 'A2': 0, 'E': 2}
    x = [np.trace(perm_matrix(p, 2)) for p in ({0: 0, 1: 1}, {0: 1, 1: 0})]
    sym, anti = int(x[0] + x[1]) // 2, int(x[0] - x[1]) // 2
    ok &= (sym, anti) == (3, 1)
    print('methyl:', methyl, ' amino: symmetric %d, antisymmetric %d' % (sym, anti))
    e5 = {i: i for i in range(5)}
    t5 = {0: 1, 1: 2, 2: 0, 3: 3, 4: 4}
    b5 = {0: 0, 1: 2, 2: 1, 3: 4, 4: 3}
    protons = [perm_matrix(p, 5) for p in (e5, t5, b5)]
    five = decompose([np.trace(P) for P in protons])
    ok &= five == {'A1': 12, 'A2': 4, 'E': 8}
    # the product of the fragments: b swaps 4 and 5, so the amino symmetric states go with A1 and the antisymmetric
    # ones exchange A1 and A2
    prod = {'A1': methyl['A1'] * sym + methyl['A2'] * anti, 'A2': methyl['A1'] * anti + methyl['A2'] * sym,
            'E': methyl['E'] * (sym + anti)}
    ok &= prod == five
    print('five protons:', five, ' product of the fragments:', prod)
    full = decompose([np.trace(np.kron(P, np.eye(1 * 3))) for P in protons])   # x C^1 (carbon-12) x C^3 (nitrogen-14)
    ok &= full == {'A1': 36, 'A2': 12, 'E': 24}
    weights = {}
    for parity, chi in (('even', TAB['A1']), ('odd', TAB['A2'])):
        weights[parity] = {g: sum(full[s] * inner(TAB[s], [a * c for a, c in zip(TAB[g], chi)]) for s in TAB) for g in TAB}
    ok &= weights == {'even': {'A1': 36, 'A2': 12, 'E': 24}, 'odd': {'A1': 12, 'A2': 36, 'E': 24}}
    print('full spin space:', full, ' weights:', weights)
    if not ok:
        sys.exit('verify_spin: mismatch')
    print('verify_spin: all checks pass')

if __name__ == '__main__':
    main()
