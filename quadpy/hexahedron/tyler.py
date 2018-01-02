# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import fs_r00, pm_rrr, z
from ..helpers import untangle


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y

        if index == 1:
            self.degree = 3
            data = [
                (frac(1, 6), fs_r00(1)),
                ]
        else:
            assert index == 2
            self.degree = 5
            data = [
                (-frac(62, 45), z()),
                (frac(16, 45), fs_r00(frac(1, 2))),
                (frac(1, 45), fs_r00(1)),
                (frac(1, 72), pm_rrr(1)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 8
        return
