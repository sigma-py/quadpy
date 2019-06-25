# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ._helpers import zero, symm_s, symm_r0, concat, QuadrilateralScheme
from ..helpers import article

citation = article(
    authors=["G.W. Tyler"],
    title="Numerical integration of functions of several variables",
    journal="Canad. J. Math.",
    volume="5",
    year="1953",
    pages="393-412",
    url="https://doi.org/10.4153/CJM-1953-044-1",
)


def tyler_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(
        zero(-frac(28, 45)),
        symm_s([frac(1, 36), 1]),
        symm_r0([frac(1, 45), 1], [frac(16, 45), frac(1, 2)]),
    )
    weights *= 4
    return QuadrilateralScheme("Tyler 1", weights, points, 5, citation)


def tyler_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pm = numpy.array([+1, -1])
    r = sqrt(frac(6, 7))
    s, t = sqrt((114 - pm * 3 * sqrt(583)) / 287)
    B1 = frac(49, 810)
    B2, B3 = (178981 + pm * 2769 * sqrt(583)) / 1888920
    weights, points = concat(symm_r0([B1, r]), symm_s([B2, s], [B3, t]))
    weights *= 4
    return QuadrilateralScheme("Tyler 2", weights, points, 7, citation)


def tyler_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(
        zero(frac(449, 315)),
        symm_r0(
            [frac(37, 1260), 1], [frac(3, 28), frac(2, 3)], [-frac(69, 140), frac(1, 3)]
        ),
        symm_s([frac(7, 540), 1], [frac(32, 135), frac(1, 2)]),
    )
    weights *= 4
    return QuadrilateralScheme("Tyler 3", weights, points, 7, citation)
