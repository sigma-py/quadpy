# -*- coding: utf-8 -*-
#
"""
Preston C. Hammer and Arthur H. Stroud,
Numerical Evaluation of Multiple Integrals II,
Math. Comp. 12 (1958), 272-280,
<https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
"""
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, fsd, z
from .helpers import volume_unit_ball, NBallScheme


def hammer_stroud_11n(n, alpha, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = sqrt(frac(n + alpha, n + alpha + 2))
    data = [(frac(1, 2 * n), fsd(n, (r, 1)))]

    points, weights = untangle(data)
    weights *= volume_unit_ball(n, symbolic=symbolic)
    return NBallScheme("Hammer-Stroud 11n", n, 3, weights, points)


def hammer_stroud_12n(n, alpha, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = sqrt(frac(3 * (n + alpha + 2), (n + 2) * (n + alpha + 4)))
    B1 = frac(
        (4 - n) * (n + 2) * (n + alpha) * (n + alpha + 4), 18 * n * (n + alpha + 2) ** 2
    )
    B2 = frac((n + 2) * (n + alpha) * (n + alpha + 4), 36 * n * (n + alpha + 2) ** 2)
    B0 = 1 - 2 * n * B1 - 2 * n * (n - 1) * B2

    data = [(B0, z(n)), (B1, fsd(n, (r, 1))), (B2, fsd(n, (r, 2)))]

    points, weights = untangle(data)
    weights *= volume_unit_ball(n, symbolic=symbolic)
    return NBallScheme("Hammer-Stroud 12n", n, 5, weights, points)


HammerStroud = {"11n": hammer_stroud_11n, "12n": hammer_stroud_12n}
