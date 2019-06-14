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

from .helpers import BallScheme
from ..helpers import untangle, fsd, pm, z


def hammer_stroud_11_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    data = [(frac(1, 6), fsd(3, (sqrt(frac(3, 5)), 1)))]
    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    return BallScheme("Hammer-Stroud 11-3", 3, weights, points)


def hammer_stroud_12_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    alpha = sqrt(frac(3, 7))
    data = [
        (frac(1, 15), z(3)),
        (frac(7, 90), fsd(3, (alpha, 1))),
        (frac(7, 180), fsd(3, (alpha, 2))),
    ]
    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    return BallScheme("Hammer-Stroud 12-3", 5, weights, points)


def hammer_stroud_14_3(variant_a=True, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    t = 1 if variant_a else -1

    sqrt14 = sqrt(14)

    # ERR The article falsely gives 0.50824... instead of 0.050824...
    a1 = frac(1, 125) * (9 + t * 2 * sqrt14)
    c1 = (71 - t * 12 * sqrt14) / 1000

    nu = sqrt((7 - t * sqrt14) / 7)
    eta1 = sqrt(5 / (21 - t * 2 * sqrt14))

    data = [(a1, fsd(3, (nu, 1))), (c1, pm(3, eta1))]

    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    name = "Hammer-Stroud 14-3" + ("a" if variant_a else "b")
    return BallScheme(name, 5, weights, points)


def hammer_stroud_15_3(variant_a=True, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    t = 1 if variant_a else -1

    sqrt30 = sqrt(30)
    nu2 = (45 - t * sqrt30) / 57
    xi2 = (18 + t * sqrt30) / 42
    eta2 = 7 / (27 + t * 2 * sqrt30)

    # The extract expressions are from Stroud's book.
    a1 = 1 / nu2 ** 3 / 63
    b1 = 1 / xi2 ** 3 / 630
    c1 = 1 / eta2 ** 3 / 2520
    a0 = 1 - 6 * a1 - 12 * b1 - 8 * c1

    data = [
        (a0, z(3)),
        (a1, fsd(3, (sqrt(nu2), 1))),
        (b1, fsd(3, (sqrt(xi2), 2))),
        (c1, pm(3, sqrt(eta2))),
    ]
    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    name = "Hammer-Stroud 15-3" + ("a" if variant_a else "b")
    return BallScheme(name, 7, weights, points)


HammerStroud = {
    "11-3": hammer_stroud_11_3,
    "12-3": hammer_stroud_12_3,
    "14-3": hammer_stroud_14_3,
    "15-3": hammer_stroud_15_3,
}
