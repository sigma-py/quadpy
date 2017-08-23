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
        else:
            assert degree == 3

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

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
