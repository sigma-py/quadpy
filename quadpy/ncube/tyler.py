# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from ..helpers import untangle, fsd, z


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self, n):
        self.degree = 3

        data = [
            (fr(3-n, 3), z(n)),
            (fr(1, 6), fsd(n, (1, 1))),
            ]

        self.points, self.weights = untangle(data)
        reference_volume = 2**n
        self.weights *= reference_volume
        return
