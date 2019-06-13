# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import concat, symm_r0, pm2


class Phillips(object):
    """
    G.M. Phillips,
    Numerical integration in two and three dimensions,
    Comput J (1967) 10 (2): 202-204,
    <https://doi.org/10.1093/comjnl/10.2.202>.

    Abstract:
    Gaussian-type quadrature formulae are derived for a rectangular region of
    two or three dimensions.
    """

    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        pm = numpy.array([+1, -1])

        self.name = "Phillips"

        c = 3 * sqrt(385)
        r, s = sqrt((105 + pm * c) / 140)
        t = sqrt(frac(3, 5))

        B1, B2 = (77 - pm * c) / 891
        B3 = frac(25, 324)

        self.degree = 7
        self.weights, self.points = concat(symm_r0([B1, r], [B2, s]), pm2([B3, t, t]))
        self.weights *= 4
        return
