# -*- coding: utf-8 -*-
#
from __future__ import division

import warnings

import numpy
import sympy

from ._rabinowitz_richter import (
    rabinowitz_richter_1 as stroud_9_1,
    rabinowitz_richter_2 as stroud_11_1,
    rabinowitz_richter_3 as stroud_11_2,
    rabinowitz_richter_4 as stroud_13_1,
    rabinowitz_richter_5 as stroud_15_1,
)
from ._stroud_secrest import (
    stroud_secrest_v as stroud_5_1,
    stroud_secrest_vi as stroud_7_1,
)

from ._helpers import E2r2Scheme
from ..helpers import untangle, fsd, pm, book

_citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


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
    return E2r2Scheme("Stroud 4-1", weights, points, 4, _citation)


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
    return E2r2Scheme("Stroud 5-2", weights, points, 5, _citation)


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
    return E2r2Scheme("Stroud 7-2", weights, points, 1, _citation)


__all__ = [
    "stroud_4_1",
    "stroud_5_1",
    "stroud_5_2",
    "stroud_7_1",
    "stroud_7_2",
    "stroud_9_1",
    "stroud_11_1",
    "stroud_11_2",
    "stroud_13_1",
    "stroud_15_1",
]
