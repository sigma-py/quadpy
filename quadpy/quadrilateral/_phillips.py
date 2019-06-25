# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ..helpers import article
from ._helpers import QuadrilateralScheme, concat, pm2, symm_r0

citation = article(
    authors=["G.M. Phillips"],
    title="Numerical integration in two and three dimensions",
    journal="Comput J",
    year="1967",
    volume="10",
    number="2",
    pages="202-204",
    url="https://doi.org/10.1093/comjnl/10.2.202",
)
# TODO add scheme for hex


def phillips(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pm = numpy.array([+1, -1])

    c = 3 * sqrt(385)
    r, s = sqrt((105 + pm * c) / 140)
    t = sqrt(frac(3, 5))

    B1, B2 = (77 - pm * c) / 891
    B3 = frac(25, 324)

    weights, points = concat(symm_r0([B1, r], [B2, s]), pm2([B3, t, t]))
    weights *= 4
    return QuadrilateralScheme("Phillips", weights, points, 7, citation)
