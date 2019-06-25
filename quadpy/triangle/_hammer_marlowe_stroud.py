# -*- coding: utf-8 -*-
#
"""
Two of the schemes also appear in

P.C. Hammer, Arthur H. Stroud,
Numerical Evaluation of Multiple Integrals II,
Mathematical Tables and Other Aids to Computation.
Vol. 12, No. 64 (Oct., 1958), pp. 272-280,
<https://www.jstor.org/stable/2002370>
"""
from __future__ import division

import numpy
import sympy

from ._helpers import TriangleScheme, s3, r, concat
from ..helpers import article

citation = article(
    authors=["P.C. Hammer", "O.J. Marlowe", "A.H. Stroud"],
    title="Numerical Integration Over Simplexes and Cones",
    journal="Mathematical Tables and Other Aids to Computation",
    volume="10",
    number="55",
    month="jul",
    year="1956",
    pages="130-137",
    url="https://doi.org/10.1090/S0025-5718-1956-0086389-6",
)


def hammer_marlowe_stroud_1(symbolic=False):
    weights, points = s3(1)
    return TriangleScheme("Hammer-Marlowe-Stroud 1", weights, points, 1, citation)


def hammer_marlowe_stroud_2(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, points = r([frac(1, 3), frac(1, 2)])
    return TriangleScheme("Hammer-Marlowe-Stroud 2", weights, points, 2, citation)


def hammer_marlowe_stroud_3(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, points = r([frac(1, 3), -frac(1, 2)])
    return TriangleScheme("Hammer-Marlowe-Stroud 3", weights, points, 2, citation)


def hammer_marlowe_stroud_4(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, points = concat(s3(-frac(9, 16)), r([frac(25, 48), frac(2, 5)]))
    return TriangleScheme("Hammer-Marlowe-Stroud 4", weights, points, 3, citation)


def hammer_marlowe_stroud_5(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    w1, w2 = [(155 - i * sqrt(15)) / 1200 for i in [+1, -1]]
    x1, x2 = [(1 + i * sqrt(15)) / 7 for i in [+1, -1]]
    weights, points = concat(s3(frac(9, 40)), r([w1, x1], [w2, x2]))
    return TriangleScheme("Hammer-Marlowe-Stroud 5", weights, points, 5, citation)
