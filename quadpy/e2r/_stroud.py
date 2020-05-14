import numpy
from sympy import Rational as frac
from sympy import cos, pi, sin, sqrt

from ..helpers import book, untangle
from ._helpers import E2rScheme
from ._rabinowitz_richter import (
    rabinowitz_richter_1 as stroud_9_1,
)  # ERR misprint in Stroud copied from original article; rabinowitz_richter_4 as stroud_13_1,
from ._rabinowitz_richter import rabinowitz_richter_2 as stroud_11_1
from ._rabinowitz_richter import rabinowitz_richter_3 as stroud_11_2
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
    return E2rScheme("Stroud 4-1", weights, points, 4, _source)


__all__ = [
    "stroud_4_1",
    "stroud_5_1",
    "stroud_7_1",
    "stroud_9_1",
    "stroud_11_1",
    "stroud_11_2",
    "stroud_15_1",
]
