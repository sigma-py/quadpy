# -*- coding: utf-8 -*-
#
"""
J. Albrecht, L. Collatz,
Zur numerischen Auswertung mehrdimensionaler Integrale,
ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
<https://doi.org/10.1002/zamm.19580380102>
"""
from __future__ import division

import numpy
import sympy

from .helpers import concat, symm_r0, symm_s, pm2, pm, zero, QuadrilateralScheme


def albrecht_collatz_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    name = "Albrecht-Collatz 1"
    degree = 3
    weights, points = concat(
        zero(frac(5, 12)), symm_r0([frac(1, 8), 1]), symm_s([frac(1, 48), 1])
    )
    weights *= 4
    return QuadrilateralScheme(name, degree, weights, points)


def albrecht_collatz_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    name = "Albrecht-Collatz 2"
    degree = 5
    r = sqrt(frac(3, 5))
    s = sqrt(frac(1, 3))
    t = sqrt(frac(14, 15))
    weights, points = concat(
        zero(frac(2, 7)), pm([frac(5, 63), 0, t]), pm2([frac(5, 36), r, s])
    )
    weights *= 4
    return QuadrilateralScheme(name, degree, weights, points)


def albrecht_collatz_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    name = "Albrecht-Collatz 3"
    degree = 5

    r = sqrt(frac(7, 15))
    s, t = [sqrt((7 + i * sqrt(24)) / 15) for i in [+1, -1]]
    weights, points = concat(
        zero(frac(2, 7)),
        pm([frac(25, 168), r, r], [frac(5, 48), +s, -t], [frac(5, 48), +t, -s]),
    )
    weights *= 4
    return QuadrilateralScheme(name, degree, weights, points)


def albrecht_collatz_4(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    name = "Albrecht-Collatz 4"
    degree = 5

    weights, points = concat(
        zero(frac(2, 45)),
        symm_r0([frac(2, 45), 1]),
        symm_s([frac(1, 60), 1], [frac(8, 45), frac(1, 2)]),
    )
    weights *= 4
    return QuadrilateralScheme(name, degree, weights, points)
