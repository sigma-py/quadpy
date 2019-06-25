# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ._helpers import BallScheme
from ..helpers import untangle, fsd, pm, z, article

_citation = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    number="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)


def hammer_stroud_11_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    data = [(frac(1, 6), fsd(3, (sqrt(frac(3, 5)), 1)))]
    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    return BallScheme("Hammer-Stroud 11-3", _citation, 3, weights, points)


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
    return BallScheme("Hammer-Stroud 12-3", _citation, 5, weights, points)


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
    return BallScheme(name, _citation, 5, weights, points)


def _hammer_stroud_15_3(variant_a, symbolic=False):
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
    return BallScheme(name, _citation, 7, weights, points)


def hammer_stroud_15_3a(symbolic=False):
    return _hammer_stroud_15_3(True, symbolic)


def hammer_stroud_15_3b(symbolic=False):
    return _hammer_stroud_15_3(False, symbolic)
