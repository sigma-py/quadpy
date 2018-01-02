# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _symm_r_0, _pm2
from ..helpers import untangle


class Phillips(object):
    '''
    G.M. Phillips,
    Numerical integration in two and three dimensions,
    Comput J (1967) 10 (2): 202-204,
    <https://doi.org/10.1093/comjnl/10.2.202>.

    Abtract:
    Gaussian-type quadrature formulae are derived for a rectangular region of
    two or three dimensions.
    '''
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        pm = numpy.array([+1, -1])

        self.name = 'Phillips'

        c = 3*sqrt(385)
        r, s = sqrt((105 + pm*c) / 140)
        t = sqrt(frac(3, 5))

        B1, B2 = (77 - pm*c) / 891
        B3 = frac(25, 324)

        self.degree = 7
        data = [
            (B1, _symm_r_0(r)),
            (B2, _symm_r_0(s)),
            (B3, _pm2(t, t))
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
