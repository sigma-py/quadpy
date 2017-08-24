# -*- coding: utf-8 -*-
#
from __future__ import division

from ..helpers import untangle, rd


class Lauffer(object):
    '''
    R. Lauffer,
    Interpolation mehfacher Integrale,
    Arch. Math. v. 6. 1955, pp. 159-164, MR 16, 862,
    <https://doi.org/10.1007/BF01900222>.
    '''
    def __init__(self, n, degree):
        self.dim = n
        self.degree = degree
        if degree == 1:
            data = [
                (1.0/(n+1), rd(n+1, [(1.0, 1)]))
                ]
        elif degree == 2:
            B = (2-n) / (n+1) / (n+2)
            C = 4 / (n+1) / (n+2)

            data = [
                (B, rd(n+1, [(1.0, 1)])),
                (C, rd(n+1, [(0.5, 2)])),
                ]
        elif degree == 3:

            B = (n**2-4*n+6) / (n+1) / (n+2) / (n+3)
            C = (27-9*n) / 2 / (n+1) / (n+2) / (n+3)
            D = 27 / (n+1) / (n+2) / (n+3)

            r = 1/3
            s = 2/3

            data = [
                (B, rd(n+1, [(1.0, 1)])),
                (C, rd(n+1, [(r, 1), (s, 1)])),
                (D, rd(n+1, [(r, 3)])),
                ]
        elif degree == 4:
            assert n >= 3

            nprod = (n+1) * (n+2) * (n+3) * (n+4)
            B1 = (-3*n**3 + 17*n**2 - 58*n + 72) / 3 / nprod
            B2 = 16 * (n**2 - 5*n + 12) / 3 / nprod
            B3 = 4 * (n**2 - 9*n + 12) / nprod
            B4 = 64 * (4-n) / 2 / nprod
            B5 = 256 / nprod

            r = 1/4
            s = 3/4
            t = 1/2

            data = [
                (B1, rd(n+1, [(1.0, 1)])),
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
            B1 = (12*n**4 - 82*n**3 + 477*n**2 - 1277*n + 1440) / 12 / nprod
            B2 = 5**2 * (-3*n**3 + 19*n**2 - 96*n + 170) / 12 / nprod
            B3 = 5**2 * (-n**3 + 13*n**2 - 47*n + 65) / 6 / nprod
            B4 = 5**3 * (n**2 - 6*n + 20) / 3 / nprod
            B5 = 5**3 * (n**2 - 11*n + 20) / 4 / nprod
            B6 = 5**4 * (5 - n) / 2 / nprod
            B7 = 5**5 / nprod

            r = 1/5
            s = 4/5
            u = 2/5
            v = 3/5

            data = [
                (B1, rd(n+1, [(1.0, 1)])),
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
