# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from ..helpers import untangle, rd


class Lauffer(object):
    '''
    R. Lauffer,
    Interpolation mehfacher Integrale,
    Arch. Math. v. 6. 1955, pp. 159-164, MR 16, 862,
    <https://doi.org/10.1007/BF01900222>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, degree):
        self.dim = n
        self.degree = degree
        if degree == 1:
            data = [
                (fr(1, n+1), rd(n+1, [(1, 1)]))
                ]
        elif degree == 2:
            B = fr(2-n, (n+1) * (n+2))
            C = fr(4, (n+1) * (n+2))

            data = [
                (B, rd(n+1, [(1, 1)])),
                (C, rd(n+1, [(fr(1, 2), 2)])),
                ]
        elif degree == 3:

            B = fr(n**2-4*n+6, (n+1) * (n+2) * (n+3))
            C = fr(27-9*n, 2 * (n+1) * (n+2) * (n+3))
            D = fr(27, (n+1) * (n+2) * (n+3))

            r = fr(1, 3)
            s = fr(2, 3)

            data = [
                (B, rd(n+1, [(1, 1)])),
                (C, rd(n+1, [(r, 1), (s, 1)])),
                (D, rd(n+1, [(r, 3)])),
                ]
        elif degree == 4:
            assert n >= 3

            nprod = (n+1) * (n+2) * (n+3) * (n+4)
            B1 = fr(-3*n**3 + 17*n**2 - 58*n + 72, 3*nprod)
            B2 = fr(16 * (n**2 - 5*n + 12), 3*nprod)
            B3 = fr(4 * (n**2 - 9*n + 12), nprod)
            B4 = fr(64 * (4-n), 2*nprod)
            B5 = fr(256, nprod)

            r = fr(1, 4)
            s = fr(3, 4)
            t = fr(1, 2)

            data = [
                (B1, rd(n+1, [(1, 1)])),
                (B2, rd(n+1, [(r, 1), (s, 1)])),
                (B3, rd(n+1, [(t, 2)])),
                (B4, rd(n+1, [(r, 2), (t, 1)])),
                (B5, rd(n+1, [(r, 4)])),
                ]

        else:
            assert degree == 5
            assert n >= 4

            nprod = (n+1) * (n+2) * (n+3) * (n+4) * (n+5)

            # ERR Stroud is missing the factor 1/12 in B1.
            B1 = fr((12*n**4 - 82*n**3 + 477*n**2 - 1277*n + 1440), 12*nprod)
            B2 = fr(5**2 * (-3*n**3 + 19*n**2 - 96*n + 170), 12*nprod)
            B3 = fr(5**2 * (-n**3 + 13*n**2 - 47*n + 65), 6*nprod)
            B4 = fr(5**3 * (n**2 - 6*n + 20), 3*nprod)
            B5 = fr(5**3 * (n**2 - 11*n + 20), 4*nprod)
            B6 = fr(5**4 * (5 - n), 2*nprod)
            B7 = fr(5**5, nprod)

            r = fr(1, 5)
            s = fr(4, 5)
            u = fr(2, 5)
            v = fr(3, 5)

            data = [
                (B1, rd(n+1, [(1, 1)])),
                (B2, rd(n+1, [(r, 1), (s, 1)])),
                (B3, rd(n+1, [(u, 1), (v, 1)])),
                (B4, rd(n+1, [(r, 2), (v, 1)])),
                (B5, rd(n+1, [(r, 1), (u, 2)])),
                (B6, rd(n+1, [(r, 3), (u, 1)])),
                (B7, rd(n+1, [(r, 5)])),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
