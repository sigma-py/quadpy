# -*- coding: utf-8 -*-
#
"""
A.H. Stroud and D. Secrest,
Approximate integration formulas for certain spherically symmetric regions,
Math. Comp. 17 (1963), 105-135,
<https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
"""
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, pm_array, pm, fsd


def v(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    nu = 2 * sqrt(5)
    xi = sqrt(5)
    eta = sqrt(15)

    data = [
        (frac(7, 10), numpy.array([[0, 0]])),
        (frac(1, 20), numpy.array([[+nu, 0], [-nu, 0]])),
        (frac(1, 20), pm_array([xi, eta])),
    ]
    return 5, data


def vi(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt74255 = sqrt(74255)

    nu = sqrt(42)
    xi, eta = [sqrt((6615 - p_m * 21 * sqrt74255) / 454) for p_m in [+1, -1]]
    A = frac(5, 588)
    B, C = [(5272105 + p_m * 18733 * sqrt74255) / 43661940 for p_m in [+1, -1]]

    data = [(A, fsd(2, (nu, 1))), (B, pm(2, xi)), (C, pm(2, eta))]
    return 7, data


_gen = {"V": v, "VI": vi}


class StroudSecrest(object):
    keys = _gen.keys()

    def __init__(self, key, symbolic=False):
        self.degree, data = _gen[key](symbolic)
        self.points, self.weights = untangle(data)
        pi = sympy.pi if symbolic else numpy.pi
        self.weights *= 2 * pi
        return
