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
    def __init__(self, n, index):
        self.dim = n
        if index == 1:
            self.degree = 1
            data = [
                (1.0/(n+1), rd(n+1, [(1.0, 1)]))
                ]
        else:
            assert index == 2
            self.degree = 2

            B = (2-n) / (n+1) / (n+2)
            C = 4 / (n+1) / (n+2)

            data = [
                (B, rd(n+1, [(1.0, 1)])),
                (C, rd(n+1, [(0.5, 2)])),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
