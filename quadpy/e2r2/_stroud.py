import warnings

import numpy
from sympy import Rational as frac
from sympy import cos, sin, sqrt

from ..helpers import book, fsd, pm, untangle
from ._helpers import E2r2Scheme
from ._rabinowitz_richter import rabinowitz_richter_1 as stroud_9_1
from ._rabinowitz_richter import rabinowitz_richter_2 as stroud_11_1
from ._rabinowitz_richter import rabinowitz_richter_3 as stroud_11_2
from ._rabinowitz_richter import rabinowitz_richter_4 as stroud_13_1
from ._rabinowitz_richter import rabinowitz_richter_5 as stroud_15_1
from ._stroud_secrest import stroud_secrest_5 as stroud_5_1
from ._stroud_secrest import stroud_secrest_6 as stroud_7_1

_source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_4_1():
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
    return E2r2Scheme("Stroud 4-1", weights, points, 4, _source)


def stroud_5_2():
    # Cartesian product Gauss formula
    r = sqrt(frac(3, 2))
    data = [
        (frac(4, 9), [[0, 0]]),
        (frac(1, 9), fsd(2, (r, 1))),
        (frac(1, 36), pm([r, r])),
    ]

    points, weights = untangle(data)
    return E2r2Scheme("Stroud 5-2", weights, points, 5, _source)


def stroud_7_2():
    # Cartesian product Gauss formula
    sqrt6 = sqrt(6)
    r, s = [sqrt((3 + p_m * sqrt6) / 2) for p_m in [+1, -1]]
    A, B = [(5 - p_m * 2 * sqrt6) / 48 for p_m in [+1, -1]]
    C = frac(1, 48)

    data = [(A, fsd(2, (r, 1))), (B, fsd(2, (s, 1))), (C, fsd(2, (r, 1), (s, 1)))]

    points, weights = untangle(data)

    # TODO find what's wrong
    warnings.warn("Stroud's Gauss product formula has degree 1, not 7.")
    return E2r2Scheme("Stroud 7-2", weights, points, 1, _source)


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
