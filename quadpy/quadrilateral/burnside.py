# -*- coding: utf-8 -*-
#
"""
W. Burnside,
An approximate quadrature formula,
Messenger of Math., v. 37, 1908, pp. 166-167.
"""
from __future__ import division

import numpy
import sympy

from .helpers import concat, symm_r0, symm_s, QuadrilateralScheme


def Burnside(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    name = "Burnside"
    degree = 5
    r = sqrt(frac(7, 15))
    s = sqrt(frac(7, 9))

    weights, points = concat(symm_r0([frac(10, 49), r]), symm_s([frac(9, 196), s]))
    weights *= 4
    return QuadrilateralScheme(name, degree, weights, points)
