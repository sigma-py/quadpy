# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ..helpers import book, untangle
from ._helpers import E2rScheme
from ._rabinowitz_richter import (
    rabinowitz_richter_1 as stroud_9_1,
)  # ERR misprint in Stroud copied from original article; rabinowitz_richter_4 as stroud_13_1,
from ._rabinowitz_richter import rabinowitz_richter_2 as stroud_11_1
from ._rabinowitz_richter import rabinowitz_richter_3 as stroud_11_2
from ._rabinowitz_richter import rabinowitz_richter_5 as stroud_15_1
from ._stroud_secrest import stroud_secrest_v as stroud_5_1
from ._stroud_secrest import stroud_secrest_vi as stroud_7_1

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
        2
        * sqrt(5)
        * numpy.array(
            [
                [cos(2 * i * pi / 5) for i in range(5)],
                [sin(2 * i * pi / 5) for i in range(5)],
            ]
        ).T
    )
    data = [(frac(7, 10), numpy.array([[0, 0]])), (frac(3, 50), pts)]

    points, weights = untangle(data)
    weights *= 2 * pi
    return E2rScheme("Stroud 4-1", weights, points, 4, _citation)


__all__ = [
    "stroud_4_1",
    "stroud_5_1",
    "stroud_7_1",
    "stroud_9_1",
    "stroud_11_1",
    "stroud_11_2",
    "stroud_15_1",
]
