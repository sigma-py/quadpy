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

from .helpers import E3r2Scheme
from ..helpers import untangle, pm_roll, fsd, pm


def stroud_secrest_vii(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    # article:
    # nu, xi = numpy.sqrt((15 + plus_minus * 3*numpy.sqrt(5)))
    # A = 3/5
    # B = 1/30

    # book:
    nu, xi = [sqrt((5 - p_m * sqrt(5)) / 4) for p_m in [+1, -1]]
    A = frac(2, 5)
    B = frac(1, 20)

    data = [(A, numpy.array([[0, 0, 0]])), (B, pm_roll(3, [nu, xi]))]
    points, weights = untangle(data)
    weights *= sqrt(pi) ** 3
    return E3r2Scheme("Stroud-Secrest VII", 5, weights, points)


def stroud_secrest_viii_a(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    r = sqrt(frac(5, 4))
    s = sqrt(frac(5, 2))
    data = [(frac(4, 25), fsd(3, (r, 1))), (frac(1, 200), pm(3, s))]
    points, weights = untangle(data)
    weights *= sqrt(pi) ** 3
    return E3r2Scheme("Stroud-Secrest VIIIa", 5, weights, points)


def stroud_secrest_viii_b(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    r = sqrt(frac(5, 2))
    s = sqrt(frac(5, 6))
    data = [
        (frac(2, 5), numpy.array([[0, 0, 0]])),
        (frac(1, 25), fsd(3, (r, 1))),
        (frac(9, 200), pm(3, s)),
    ]
    points, weights = untangle(data)
    weights *= sqrt(pi) ** 3
    return E3r2Scheme("Stroud-Secrest VIIIb", 5, weights, points)


def stroud_secrest_ix(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    r, s = [sqrt((15 - p_m * 5 * sqrt(5)) / 12) for p_m in [+1, -1]]
    t = sqrt(frac(5, 6))

    data = [
        (frac(2, 5), numpy.array([[0, 0, 0]])),
        (frac(3, 100), pm_roll(3, [r, s])),
        (frac(3, 100), pm(3, t)),
    ]
    points, weights = untangle(data)
    weights *= sqrt(pi) ** 3
    return E3r2Scheme("Stroud-Secrest IX", 5, weights, points)


def stroud_secrest_x(positive, symbolic=False):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    plus_minus = 1 if positive else -1

    sqrt15 = sqrt(15)

    r = sqrt((15 + plus_minus * sqrt15) / 4)
    s = sqrt((6 - plus_minus * sqrt15) / 2)
    t = sqrt((9 + plus_minus * 2 * sqrt15) / 2)
    A = (720 + plus_minus * 8 * sqrt15) / 2205
    B = (270 - plus_minus * 46 * sqrt15) / 15435
    C = (162 + plus_minus * 41 * sqrt15) / 6174
    D = (783 - plus_minus * 202 * sqrt15) / 24696

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, fsd(3, (r, 1))),
        (C, fsd(3, (s, 2))),
        (D, pm(3, t)),
    ]
    points, weights = untangle(data)
    weights *= sqrt(pi) ** 3
    return E3r2Scheme("Stroud-Secrest X", 7, weights, points)


def stroud_secrest_xi(positive, symbolic=False):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    p_m = 1 if positive else -1

    sqrt2 = sqrt(2)
    sqrt5 = sqrt(5)
    sqrt10 = sqrt(10)

    r = sqrt((25 + p_m * 15 * sqrt2 + 5 * sqrt5 + p_m * 3 * sqrt10) / 4)
    s = sqrt((25 + p_m * 15 * sqrt2 - 5 * sqrt5 - p_m * 3 * sqrt10) / 4)
    t = sqrt((3 - p_m * sqrt2) / 2)
    u = sqrt((9 - p_m * 3 * sqrt2 - 3 * sqrt5 + p_m * sqrt10) / 4)
    v = sqrt((9 - p_m * 3 * sqrt2 + 3 * sqrt5 - p_m * sqrt10) / 4)

    A = (80 + p_m * 8 * sqrt2) / 245
    B = (395 - p_m * 279 * sqrt2) / 13720
    C = (45 + p_m * 29 * sqrt2) / 2744

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, pm_roll(3, [r, s])),
        (C, pm_roll(3, [u, v])),
        (C, pm(3, t)),
    ]
    points, weights = untangle(data)
    weights *= sqrt(pi) ** 3
    return E3r2Scheme("Stroud-Secrest X", 7, weights, points)


StroudSecrest = {
    "VII": stroud_secrest_vii,
    "VIIIa": stroud_secrest_viii_a,
    "VIIIb": stroud_secrest_viii_b,
    "IX": stroud_secrest_ix,
    "Xa": lambda symbolic=False: stroud_secrest_x(True, symbolic),
    "Xb": lambda symbolic=False: stroud_secrest_x(False, symbolic),
    "XIa": lambda symbolic=False: stroud_secrest_xi(True, symbolic),
    "XIb": lambda symbolic=False: stroud_secrest_xi(False, symbolic),
}
