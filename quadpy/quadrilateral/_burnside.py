# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._helpers import concat, symm_r0, symm_s, QuadrilateralScheme
from ..helpers import article

citation = article(
    authors=["W. Burnside"],
    title="An approximate quadrature formula",
    journal="Messenger of Math.",
    volume="37",
    year="1908",
    pages="166-167",
)


def burnside(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = sqrt(frac(7, 15))
    s = sqrt(frac(7, 9))

    weights, points = concat(symm_r0([frac(10, 49), r]), symm_s([frac(9, 196), s]))
    weights *= 4
    return QuadrilateralScheme("Burnside", weights, points, 5, citation)
