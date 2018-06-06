# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import untangle2


class HammerStroud(object):
    """
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    """

    def __init__(self, degree, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        self.name = "HammerStroud({})".format(degree)
        self.degree = degree

        data = {
            2: {"s31": [[frac(1, 4), (5 - sqrt(5)) / 20]]},
            3: {"s4": [[-frac(4, 5)]], "s31": [[+frac(9, 20), frac(1, 6)]]},
        }[degree]

        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
