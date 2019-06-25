# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._helpers import untangle2, TetrahedronScheme
from ..helpers import article

citation = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    volume="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)


def hammer_stroud_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 2
    data = {"s31": [[frac(1, 4), (5 - sqrt(5)) / 20]]}
    points, weights = untangle2(data)
    return TetrahedronScheme("Hammer-Stroud 2", weights, points, degree, citation)


def hammer_stroud_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 3
    data = {"s4": [[-frac(4, 5)]], "s31": [[+frac(9, 20), frac(1, 6)]]}
    points, weights = untangle2(data)
    return TetrahedronScheme("Hammer-Stroud 2", weights, points, degree, citation)
