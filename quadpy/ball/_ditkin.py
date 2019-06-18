# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, z, pm_array0, pm, article
from ._helpers import BallScheme

_citation = article(
    authors=["V.A. Ditkin"],
    title="On certain approximate formulas for the calculation of triple integrals",
    journal="Doklady Akad. Nauk SSSR (N.S.)",
    number="62",
    year=1948,
    pages="445–447",
    comment="Russian",
)


def ditkin_1(alpha=0, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    pm_ = numpy.array([+1, -1])

    B0 = frac(4, (alpha + 5) ** 2)
    B1 = frac((alpha + 3) * (alpha + 7), 12 * (alpha + 5) ** 2)

    r_s = sqrt((alpha + 5) * (5 + pm_ * sqrt(5)) / 10 / (alpha + 7))

    data = [
        (B0, z(3)),
        (B1, pm_array0(3, r_s, [0, 1])),
        (B1, pm_array0(3, r_s, [1, 2])),
        (B1, pm_array0(3, r_s, [2, 0])),
    ]

    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    return BallScheme("Ditkin 1", _citation, 5, weights, points)


def ditkin_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    pm_ = numpy.array([+1, -1])

    B0 = frac(4, 25)
    B1 = frac(21, 500)

    r_s = sqrt((15 + pm_ * 5 * sqrt(5)) / 42)
    t = sqrt(frac(5, 21))

    data = [
        (B0, z(3)),
        (B1, pm_array0(3, r_s, [0, 1])),
        (B1, pm_array0(3, r_s, [1, 2])),
        (B1, pm_array0(3, r_s, [2, 0])),
        (B1, pm(3, t)),
    ]

    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    return BallScheme("Ditkin 2", _citation, 5, weights, points)


def ditkin_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    pm_ = numpy.array([+1, -1])

    B0 = frac(16, 175)
    B1 = frac(81, 1400)
    B2 = frac(3, 280)

    sqrt5 = sqrt(5)
    r_s = sqrt((5 + pm_ * sqrt5) / 18)
    t = sqrt(frac(1, 3))
    u_v = sqrt((3 - pm_ * sqrt5) / 6)

    data = [
        (B0, z(3)),
        (B1, pm_array0(3, r_s, [0, 1])),
        (B1, pm_array0(3, r_s, [1, 2])),
        (B1, pm_array0(3, r_s, [2, 0])),
        (B2, pm_array0(3, u_v, [0, 1])),
        (B2, pm_array0(3, u_v, [1, 2])),
        (B2, pm_array0(3, u_v, [2, 0])),
        (B2, pm(3, t)),
    ]

    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    return BallScheme("Ditkin 3", _citation, 7, weights, points)
