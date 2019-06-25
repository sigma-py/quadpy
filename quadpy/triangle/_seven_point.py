# -*- coding: utf-8 -*-
#

import sympy

from ._helpers import concat, s2, s3, TriangleScheme


def seven_point(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = concat(
        s3(frac(9, 20)), s2([frac(1, 20), 0], [frac(2, 15), frac(1, 2)])
    )
    return TriangleScheme("Seven-point scheme", weights, points, 3)
