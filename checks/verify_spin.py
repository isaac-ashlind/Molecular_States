"""Spin-weight structure behind figure 8: the five proton spins as methyl (3) times amino (2).

Checks, by explicit permutation matrices on the 2^n spin basis:
  methyl (C^2)^3 under S3 = <(123),(23)>  decomposes as 4 A1 + 2 E   (quartet S=3/2, two doublets S=1/2)
  amino  (C^2)^2 under the swap (45)      decomposes as 3 (+) + 1 (-) (triplet, singlet)
  all five under G6 = <t=(123), (23)(45)> decomposes as 12 A1 + 4 A2 + 8 E, i.e. the product of the two.
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

CLS = [1, 2, 3]  # class sizes of S3: e, 3-cycles, transpositions
TAB = {'A1': [1, 1, 1], 'A2': [1, 1, -1], 'E': [2, -1, 0]}

def decompose(chars):
    return {k: round(sum(c * x * y for c, x, y in zip(CLS, chars, v)) / 6) for k, v in TAB.items()}

def main():
    ok = True
    e3, c3, s3 = {0: 0, 1: 1, 2: 2}, {0: 1, 1: 2, 2: 0}, {0: 0, 1: 2, 2: 1}
    methyl = decompose([np.trace(perm_matrix(p, 3)) for p in (e3, c3, s3)])
    ok &= methyl == {'A1': 4, 'A2': 0, 'E': 2}
    print('methyl 8 states:', methyl)
    e2, s2 = {0: 0, 1: 1}, {0: 1, 1: 0}
    x = [np.trace(perm_matrix(p, 2)) for p in (e2, s2)]
    sym, anti = (x[0] + x[1]) / 2, (x[0] - x[1]) / 2
    ok &= (sym, anti) == (3, 1)
    print('amino 4 states: symmetric %d, antisymmetric %d' % (sym, anti))
    e5 = {i: i for i in range(5)}
    t5 = {0: 1, 1: 2, 2: 0, 3: 3, 4: 4}
    b5 = {0: 0, 1: 2, 2: 1, 3: 4, 4: 3}
    five = decompose([np.trace(perm_matrix(p, 5)) for p in (e5, t5, b5)])
    ok &= five == {'A1': 12, 'A2': 4, 'E': 8}
    print('five protons 32 states:', five)
    prod = {'A1': methyl['A1'] * 3, 'A2': methyl['A1'] * 1, 'E': methyl['E'] * (3 + 1)}
    ok &= prod == five
    print('product of the fragments:', prod, 'matches' if prod == five else 'MISMATCH')
    # ---- the rigid limit (closing page): under H = {E, b} the spin space is 20 A' + 12 A''; each rigid species induces
    # a tunnelling multiplet of G6 (Frobenius: <Ind chi, Gamma> = <chi, Gamma|H>), A' -> A1 + E and A'' -> A2 + E, and
    # the weights add up, 20 = 12 + 8 and 12 = 4 + 8
    chi_b = np.trace(perm_matrix(b5, 5))
    rigid = {"A'": round((32 + chi_b) / 2), "A''": round((32 - chi_b) / 2)}
    ok &= rigid == {"A'": 20, "A''": 12}
    print('five protons under H = {E, b}:', rigid)
    restrict = {g: {"A'": (TAB[g][0] + TAB[g][2]) // 2, "A''": (TAB[g][0] - TAB[g][2]) // 2} for g in TAB}
    induced = {h: {g: restrict[g][h] for g in TAB} for h in ("A'", "A''")}
    ok &= induced["A'"] == {'A1': 1, 'A2': 0, 'E': 1} and induced["A''"] == {'A1': 0, 'A2': 1, 'E': 1}
    sums = {h: sum(induced[h][g] * five[g] for g in TAB) for h in induced}
    ok &= sums == rigid
    print('induced multiplets:', induced, '| weights add up:', sums)
    if not ok:
        sys.exit('verify_spin: mismatch')
    print('verify_spin: all checks pass')

if __name__ == '__main__':
    main()
