# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import unroll


class Burnside(object):
    """
    W. Burnside,
    An approximate quadrature formula,
    Messenger of Math., v. 37, 1908, pp. 166-167.
    """

    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "Burnside"
        self.degree = 5
        r = sqrt(frac(7, 15))
        s = sqrt(frac(7, 9))

        data = {
            "symm_r0": [[frac(10, 49), r]],
            "symm_s": [[frac(9, 196), s]]
        }

        self.points, self.weights = unroll(data)
        self.weights *= 4
        return
