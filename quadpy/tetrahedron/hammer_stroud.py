# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _s4, _s31
from ..helpers import untangle


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, degree, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        self.name = 'HammerStroud({})'.format(degree)
        self.degree = degree
        if degree == 2:
            data = [
                (frac(1, 4), _s31((5 - sqrt(5))/20)),
                ]
        else:
            assert degree == 3
            data = [
                (-frac(4, 5), _s4()),
                (+frac(9, 20), _s31(frac(1, 6))),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
