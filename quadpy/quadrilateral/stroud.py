# -*- coding: utf-8 -*-
#
from __future__ import division

import warnings

import numpy
import sympy

from .albrecht_collatz import (
    albrecht_collatz_1 as stroud_c2_3_4,
    albrecht_collatz_2 as stroud_c2_5_1,
    albrecht_collatz_3 as stroud_c2_5_2,
    albrecht_collatz_4 as stroud_c2_5_6,
)
from .burnside import burnside as stroud_c2_5_3
from .irwin import irwin_1 as stroud_c2_3_5, irwin_2 as stroud_c2_5_7
from .maxwell import maxwell as stroud_c2_7_3
from .meister import meister as stroud_c2_7_6
from .miller import miller as stroud_c2_1_2
from .phillips import phillips as stroud_c2_7_2
from .rabinowitz_richter import (
    rabinowitz_richter_1 as stroud_c2_9_1,
    rabinowitz_richter_2 as stroud_c2_11_1,
    rabinowitz_richter_3 as stroud_c2_11_2,
    rabinowitz_richter_4 as stroud_c2_13_1,
    rabinowitz_richter_5 as stroud_c2_15_1,
    rabinowitz_richter_6 as stroud_c2_15_2,
)
from .tyler import (
    tyler_1 as stroud_c2_5_5,
    tyler_2 as stroud_c2_7_1,
    tyler_3 as stroud_c2_7_5,
)
from .. import ncube
from .helpers import concat, zero, symm_r0, symm_s, symm_s_t, QuadrilateralScheme
from ..helpers import book

citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_c2_1_1(symbolic=False):
    # product trapezoidal
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = symm_s(frac(1, 4), 1)
    weights *= 4
    return QuadrilateralScheme("Stroud C2 1-1", weights, points, 1, citation)


def stroud_c2_3_1(symbolic=False):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    # ERR misprint in Stroud: sqrt(1/3) vs 1/3
    weights, points = symm_s([frac(1, 4), sqrt(frac(1, 3))])
    weights *= 4
    return QuadrilateralScheme("Stroud C2 3-1", weights, points, 3, citation)


def stroud_c2_3_2(symbolic=False):
    return ncube.ewing(2, symbolic)


def stroud_c2_3_3(symbolic=False):
    return ncube.stroud_cn_3_6(2, symbolic)


def stroud_c2_5_4(symbolic=False):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    r = sqrt(frac(3, 5))
    weights, points = concat(
        zero(frac(16, 81)), symm_r0([frac(10, 81), r]), symm_s([frac(25, 324), r])
    )
    weights *= 4
    return QuadrilateralScheme("Stroud C2 5-4", weights, points, 5, citation)


def stroud_c2_7_4(symbolic=False):
    # product Gauss 7
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    pm = numpy.array([+1, -1])

    r, s = sqrt((15 - pm * 2 * sqrt(30)) / 35)

    B1, B2 = (59 + pm * 6 * sqrt(30)) / 864
    B3 = frac(49, 864)

    r = sqrt(frac(3, 5))
    weights, points = concat(symm_s([B1, r], [B2, s]), symm_s_t([B3, r, s]))

    weights *= 4

    # TODO fix
    warnings.warn("Formula only has degree 1!")
    return QuadrilateralScheme("Stroud C2 7-4", weights, points, 1, citation)


__all__ = [
    "stroud_c2_1_1",
    "stroud_c2_1_2",
    "stroud_c2_3_1",
    "stroud_c2_3_2",
    "stroud_c2_3_3",
    "stroud_c2_3_4",
    "stroud_c2_3_5",
    "stroud_c2_5_1",
    "stroud_c2_5_2",
    "stroud_c2_5_3",
    "stroud_c2_5_4",
    "stroud_c2_5_5",
    "stroud_c2_5_6",
    "stroud_c2_5_7",
    "stroud_c2_7_1",
    "stroud_c2_7_2",
    "stroud_c2_7_3",
    "stroud_c2_7_4",
    "stroud_c2_7_5",
    "stroud_c2_7_6",
    "stroud_c2_9_1",
    "stroud_c2_11_1",
    "stroud_c2_11_2",
    "stroud_c2_13_1",
    "stroud_c2_15_1",
    "stroud_c2_15_2",
]
