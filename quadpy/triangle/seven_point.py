# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import concat, s2, s3, TriangleScheme


def SevenPoint(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = concat(
        s3(frac(9, 20)), s2([frac(1, 20), 0], [frac(2, 15), frac(1, 2)])
    )
    return TriangleScheme("Seven-point scheme", 3, weights, bary)
