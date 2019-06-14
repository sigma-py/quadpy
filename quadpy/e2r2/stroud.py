# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from __future__ import division

import warnings

import numpy
import sympy

from .rabinowitz_richter import RabinowitzRichter
from .stroud_secrest import StroudSecrest

from .helpers import E2r2Scheme
from ..helpers import untangle, fsd, pm


def stroud_4_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    cos = numpy.vectorize(sympy.cos) if symbolic else numpy.cos
    sin = numpy.vectorize(sympy.sin) if symbolic else numpy.sin
    pi = sympy.pi if symbolic else numpy.pi

    pts = (
        sqrt(2)
        * numpy.array(
            [
                [cos(2 * i * numpy.pi / 5) for i in range(5)],
                [sin(2 * i * numpy.pi / 5) for i in range(5)],
            ]
        ).T
    )
    data = [(frac(1, 2), numpy.array([[0, 0]])), (frac(1, 10), pts)]

    points, weights = untangle(data)
    weights *= pi
    return E2r2Scheme("Stroud 4-1", 4, weights, points)


def stroud_5_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    # Cartesian product Gauss formula
    r = sqrt(frac(3, 2))
    data = [
        (frac(4, 9), numpy.array([[0, 0]])),
        (frac(1, 9), fsd(2, (r, 1))),
        (frac(1, 36), pm(2, r)),
    ]

    points, weights = untangle(data)
    weights *= pi
    return E2r2Scheme("Stroud 5-2", 5, weights, points)


def stroud_7_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    # Cartesian product Gauss formula
    sqrt6 = sqrt(6)
    r, s = [sqrt((3 + p_m * sqrt6) / 2) for p_m in [+1, -1]]
    A, B = [(5 - p_m * 2 * sqrt6) / 48 for p_m in [+1, -1]]
    C = frac(1, 48)

    data = [(A, fsd(2, (r, 1))), (B, fsd(2, (s, 1))), (C, fsd(2, (r, 1), (s, 1)))]

    points, weights = untangle(data)
    weights *= pi

    # TODO find what's wrong
    warnings.warn("Stroud's Gauss product formula has degree 1, not 7.")
    return E2r2Scheme("Stroud 7-2", 1, weights, points)


Stroud = {
    "4-1": stroud_4_1,
    "5-1": StroudSecrest["V"],
    "5-2": stroud_5_2,
    "7-1": StroudSecrest["VI"],
    "7-2": stroud_7_2,
    "9-1": RabinowitzRichter[1],
    "11-1": RabinowitzRichter[2],
    "11-2": RabinowitzRichter[3],
    "13-1": RabinowitzRichter[4],
    "15-1": RabinowitzRichter[5],
}
