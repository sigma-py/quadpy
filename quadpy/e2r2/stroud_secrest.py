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

from .helpers import E2r2Scheme
from ..helpers import untangle, pm_array, pm_array0, fsd, pm


def stroud_secrest_v(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    nu = sqrt(2)
    xi = nu / 2
    eta = sqrt(6) / 2
    A = frac(1, 2)
    B = frac(1, 12)

    data = [
        (A, numpy.array([[0, 0]])),
        (B, pm_array0(2, [nu], [0])),
        (B, pm_array([xi, eta])),
    ]

    points, weights = untangle(data)
    weights *= pi
    return E2r2Scheme("Stroud-Secrest V", 5, weights, points)


def stroud_secrest_vi(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    sqrt5 = sqrt(5)
    nu = sqrt(3)
    xi, eta = [sqrt((9 - p_m * 3 * sqrt5) / 8) for p_m in [+1, -1]]
    A = frac(1, 36)
    B, C = [(5 + p_m * 2 * sqrt5) / 45 for p_m in [+1, -1]]

    data = [(A, fsd(2, (nu, 1))), (B, pm(2, xi)), (C, pm(2, eta))]

    points, weights = untangle(data)
    weights *= pi
    return E2r2Scheme("Stroud-Secrest VI", 7, weights, points)


StroudSecrest = {"V": stroud_secrest_v, "VI": stroud_secrest_vi}
