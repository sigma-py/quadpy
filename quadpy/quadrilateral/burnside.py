# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _symm_r_0, _symm_s
from ..helpers import untangle


class Burnside(object):
    '''
    W. Burnside,
    An approximate quadrature formula,
    Messenger of Math., v. 37, 1908, pp. 166-167.
    '''
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = 'Burnside'
        self.degree = 5
        r = sqrt(frac(7, 15))
        s = sqrt(frac(7, 9))
        data = [
            (frac(10, 49), _symm_r_0(r)),
            (frac(9, 196), _symm_s(s)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
