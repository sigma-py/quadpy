# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import fs_r00, pm_rrr, z
from ..helpers import untangle


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://dx.doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self, index):
        if index == 1:
            self.degree = 3
            data = [
                (fr(1, 6), fs_r00(1)),
                ]
        else:
            assert index == 2
            self.degree = 5
            data = [
                (-fr(62, 45), z()),
                (fr(16, 45), fs_r00(fr(1, 2))),
                (fr(1, 45), fs_r00(1)),
                (fr(1, 72), pm_rrr(1)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 8
        return
