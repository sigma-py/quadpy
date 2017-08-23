# -*- coding: utf-8 -*-
#
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
        assert index == 1
        self.degree = 1
        data = [
            (1.0/(n+1), rd(n+1, [(1.0, 1)]))
            ]
        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
