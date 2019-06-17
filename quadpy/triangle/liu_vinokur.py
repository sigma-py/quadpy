# -*- coding: utf-8 -*-
#
"""
Y. Liu and M. Vinokur,
Exact Integrations of Polynomials and Symmetric Quadrature Formulas over Arbitrary
Polyhedral Grids,
Journal of Computational Physics, 140, 122–147 (1998).
DOI: 10.1006/jcph.1998.5884,
<https://doi.org/10.1006/jcph.1998.5884>.
"""
from __future__ import division

import numpy
import sympy

from .helpers import s3, TriangleScheme, alpha, concat, gamma_delta


def liu_vinokur_1(symbolic=False):
    weights, bary = s3(1)
    return TriangleScheme("Liu-Vinokur 1", 1, weights, bary)


def liu_vinokur_2(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = alpha([frac(1, 3), 1])
    return TriangleScheme("Liu-Vinokur 2", 1, weights, bary)


def liu_vinokur_3(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = alpha([frac(1, 3), -frac(1, 2)])
    return TriangleScheme("Liu-Vinokur 3", 2, weights, bary)


def liu_vinokur_4(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = concat(s3(frac(3, 4)), alpha([frac(1, 12), 1]))
    return TriangleScheme("Liu-Vinokur 4", 2, weights, bary)


def liu_vinokur_5(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = concat(
        s3(-frac(9, 16)),
        # ERR Wrongly specified in the article as 25 (instead of 2/5).
        alpha([frac(25, 48), frac(2, 5)]),
    )
    return TriangleScheme("Liu-Vinokur 5", 3, weights, bary)


def liu_vinokur_6(symbolic=False):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt21 = sqrt(21)
    weights, bary = alpha(
        [(1 + sqrt21) / 120, 1], [(39 - sqrt21) / 120, (1 - sqrt21) / 10]
    )
    return TriangleScheme("Liu-Vinokur 6", 3, weights, bary)


def liu_vinokur_7(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = concat(
        s3(frac(9, 20)), alpha([frac(1, 20), 1], [frac(2, 15), -frac(1, 2)])
    )
    return TriangleScheme("Liu-Vinokur 7", 3, weights, bary)


def liu_vinokur_8(symbolic=False):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt10 = sqrt(10)
    sqrt_b = sqrt(950 - 220 * sqrt10)
    a1 = (-10 + 5 * sqrt10 + sqrt_b) / 30
    a2 = (-10 + 5 * sqrt10 - sqrt_b) / 30
    weights, bary = alpha(
        [(5 * a2 - 2) / (60 * a1 ** 2 * (a2 - a1)), a1],
        [(5 * a1 - 2) / (60 * a2 ** 2 * (a1 - a2)), a2],
    )
    return TriangleScheme("Liu-Vinokur 8", 4, weights, bary)


def liu_vinokur_9(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = concat(
        s3(frac(27, 80)),
        alpha([frac(8, 105), -frac(1, 2)], [frac(81, 560), frac(2, 3)]),
    )
    return TriangleScheme("Liu-Vinokur 9", 4, weights, bary)


def liu_vinokur_10(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt13 = sqrt(13)
    weights, bary = alpha(
        [(11 - 1 * sqrt13) / 360, 1],
        [(80 - 16 * sqrt13) / 360, -frac(1, 2)],
        [(29 + 17 * sqrt13) / 360, (-1 + sqrt13) / 6],
    )
    return TriangleScheme("Liu-Vinokur 10", 4, weights, bary)


def liu_vinokur_11(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt3 = sqrt(3)
    weights, bary = concat(
        s3(+frac(9, 20)),
        alpha([-frac(1, 60), 1]),
        gamma_delta([+frac(1, 10), (3 + sqrt3) / 6, (3 - sqrt3) / 6]),
    )
    return TriangleScheme("Liu-Vinokur 11", 4, weights, bary)


def liu_vinokur_12(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt15 = sqrt(15)
    weights, bary = concat(
        s3(frac(9, 40)),
        alpha(
            [(155 - sqrt15) / 1200, (1 + sqrt15) / 7],
            [(155 + sqrt15) / 1200, (1 - sqrt15) / 7],
        ),
    )
    return TriangleScheme("Liu-Vinokur 12", 5, weights, bary)


def liu_vinokur_13(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = concat(
        s3(frac(81, 320)),
        alpha(
            [frac(1, 90), 1],
            [frac(16, 225), -frac(1, 2)],
            [frac(2401, 14400), frac(4, 7)],
        ),
    )
    return TriangleScheme("Liu-Vinokur 13", 5, weights, bary)


LiuVinokur = {
    1: liu_vinokur_1,
    2: liu_vinokur_2,
    3: liu_vinokur_3,
    4: liu_vinokur_4,
    5: liu_vinokur_5,
    6: liu_vinokur_6,
    7: liu_vinokur_7,
    8: liu_vinokur_8,
    9: liu_vinokur_9,
    10: liu_vinokur_10,
    11: liu_vinokur_11,
    12: liu_vinokur_12,
    13: liu_vinokur_13,
}
