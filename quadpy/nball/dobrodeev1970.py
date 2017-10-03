# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr, gamma, pi

from ..helpers import untangle, fsd, z


class Dobrodeev1970(object):
    '''
    L.N. Dobrodeev,
    Cubature formulas of the seventh order of accuracy for a hypersphere and a
    hypercube,
    USSR Computational Mathematics and Mathematical Physics,
    Volume 10, Issue 1, 1970, Pages 252–253,
    <https://doi.org/10.1016/0041-5553(70)90084-4>.
    '''
    def __init__(self, n):
        self.name = 'Dobrodeev'
        self.degree = 7
        self.dim = n

        A = fr(1, 8)
        B = fr(5-n, 4)
        C = fr((6 - n) * (1 - n**2) + 36, 4*(n + 3))
        D = fr(81, (n + 3) * (n + 6)**2)
        E = fr(45*n**2 + 324*n + 216, n**2 + 12*n + 36) \
            - fr(n * (n**2 - 12*n + 65), 6)

        r = sqrt(fr(3, n+6))
        data = [
            (A, fsd(n, (r, 3))),
            (B, fsd(n, (r, 2))),
            (C, fsd(n, (r, 1))),
            (D, fsd(n, (1, 1))),
            (E, z(n)),
            ]

        self.points, self.weights = untangle(data)

        self.weights /= (
            fr(n, 2) * gamma(fr(n, 2)) / sqrt(pi)**n
            * fr(27 * (n+2) * (n+4), (n+6)**2)
            )
        return
