# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _symm_r_0, _symm_s, _z
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
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        pm = numpy.array([+1, -1])

        self.name = 'Tyler({})'.format(index)
        if index == 1:
            self.degree = 5
            data = [
                (-frac(28, 45), _z()),
                (frac(1, 36), _symm_s(1)),
                (frac(1, 45), _symm_r_0(1)),
                (frac(16, 45), _symm_r_0(frac(1, 2))),
                ]
        elif index == 2:
            self.degree = 7
            r = sqrt(frac(6, 7))
            s, t = sqrt((114 - pm*3*sqrt(583)) / 287)
            B1 = frac(49, 810)
            B2, B3 = (178981 + pm * 2769 * sqrt(583)) / 1888920
            data = [
                (B1, _symm_r_0(r)),
                (B2, _symm_s(s)),
                (B3, _symm_s(t)),
                ]
        else:
            assert index == 3
            self.degree = 7
            r = frac(2, 3)
            s = frac(1, 3)
            t = frac(1, 2)
            data = [
                (frac(449, 315), _z()),
                (frac(37, 1260), _symm_r_0(1)),
                (frac(3, 28), _symm_r_0(r)),
                (-frac(69, 140), _symm_r_0(s)),
                (frac(7, 540), _symm_s(1)),
                (frac(32, 135), _symm_s(t)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
