# -*- coding: utf-8 -*-
#
"""
P. Hillion,
Numerical Integration on a Triangle,
International Journal for Numerical Methods in Engineering,
Vol. 11, 797-815 (1977).
DOI:10.1002/nme.1620110504,
<https://doi.org/10.1002/nme.1620110504>.

Note that the schemes here are not fully symmetric.
"""
# ENH in the article, most schemes are given only in single precision. quadpy adds
# symbolic expressions
from __future__ import division

import numpy
import sympy

from .helpers import TriangleScheme, s3, s2, mirror, concat


def hillion_1(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = s3(frac(1, 2))
    weights *= 2
    return TriangleScheme("Hillion 1", 1, weights, bary)


def hillion_2(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = s2([frac(1, 6), frac(1, 2)])
    weights *= 2
    return TriangleScheme("Hillion 2", 2, weights, bary)


def hillion_3(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = s2([frac(1, 6), frac(1, 6)])
    weights *= 2
    return TriangleScheme("Hillion 3", 2, weights, bary)


def hillion_4(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    a0, a1 = [(3 + i * sqrt(3)) / 8 for i in [+1, -1]]
    weights, bary = concat(([frac(1, 18)], [[0, 0, 1]]), mirror([frac(2, 9), a0, a1]))
    weights *= 2
    return TriangleScheme("Hillion 4", 2, weights, bary)


def hillion_5(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    a0, a1 = [(3 + i * sqrt(3)) / 8 for i in [+1, -1]]
    weights, bary = concat(
        ([frac(1, 18)], [[frac(2, 3), frac(2, 3), -frac(1, 3)]]),
        mirror([frac(2, 9), frac(2, 3) - a0, frac(2, 3) - a1]),
    )
    weights *= 2
    return TriangleScheme("Hillion 5", 2, weights, bary)


def hillion_6(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    lm, mu = [(2 + i * sqrt(2 + i * sqrt(3))) / 6 for i in [+1, -1]]
    weights, bary = mirror(
        [frac(1, 8), lm, mu], [frac(1, 8), frac(2, 3) - lm, frac(2, 3) - mu]
    )
    weights *= 2
    return TriangleScheme("Hillion 6", 2, weights, bary)


def hillion_7(symbolic=False):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    pm = numpy.array([+1, -1])

    a, b = (6 + sqrt(2) + pm * sqrt(6 * (3 + 2 * sqrt(2)))) / 20
    c, d = (6 - sqrt(2) + pm * sqrt(6 * (3 - 2 * sqrt(2)))) / 20
    w1 = (2 - 3 * (b + c)) / 12 / (a + d - b - c)
    w2 = (2 - 3 * (a + d)) / 12 / (b + c - a - d)

    weights, bary = mirror([w1, a, d], [w2, c, b])
    weights *= 2
    return TriangleScheme("Hillion 7", 3, weights, bary)


def hillion_8(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    lambda2, lambda3 = [(32 + i * 2 * sqrt(46)) / 105 for i in [+1, -1]]
    w1, w2 = [(3266 + i * 19 * sqrt(46)) / 17664 for i in [+1, -1]]
    weights, bary = concat(
        mirror([frac(25, 384), 0, frac(4, 5)]),
        ([w1], [[lambda2, lambda2, 1 - 2 * lambda2]]),
        ([w2], [[lambda3, lambda3, 1 - 2 * lambda3]]),
    )
    weights *= 2
    return TriangleScheme("Hillion 8", 3, weights, bary)


def hillion_9(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    # ERR the article is missing the minus sign
    weights, bary = concat(s3(-frac(9, 32)), s2([frac(25, 96), frac(1, 5)]))
    weights *= 2
    return TriangleScheme("Hillion 9", 3, weights, bary)


def hillion_10(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    lambda1, lambda2 = [(16 + i * 2 * sqrt(14)) / 25 for i in [+1, -1]]
    w1, w2 = [(161 + i * 17 * sqrt(14)) / 2688 for i in [+1, -1]]
    weights, bary = concat(
        mirror([w2, lambda1, 0], [w1, 0, lambda2]),
        ([frac(25, 96)], [[frac(2, 5), frac(2, 5), frac(1, 5)]]),
    )
    weights *= 2
    return TriangleScheme("Hillion 10", 3, weights, bary)


Hillion = {
    1: hillion_1,
    2: hillion_3,
    3: hillion_4,
    4: hillion_4,
    5: hillion_5,
    6: hillion_6,
    7: hillion_7,
    8: hillion_8,
    9: hillion_9,
    10: hillion_10,
}
