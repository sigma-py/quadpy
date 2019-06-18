# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._albrecht import (
    albrecht_4 as stroud_s2_9_1,
    albrecht_5 as stroud_s2_11_2,
    albrecht_6 as stroud_s2_13_2,
    albrecht_7 as stroud_s2_15_2,
    albrecht_8 as stroud_s2_17_1,
)
from ._albrecht_collatz import albrecht_collatz as stroud_s2_3_2
from ._hammer_stroud import (
    hammer_stroud_11_2 as stroud_s2_3_1,
    hammer_stroud_12_2 as stroud_s2_5_2,
    hammer_stroud_18 as stroud_s2_7_2,
)
from ._mysovskih import (
    mysovskih_1 as stroud_s2_4_1,
    mysovskih_2 as stroud_s2_11_1,
    mysovskih_3 as stroud_s2_15_1,
)
from ._peirce_1956 import (
    peirce_1956_1 as stroud_s2_7_1,
    peirce_1956_2 as stroud_s2_9_5,
    peirce_1956_3 as stroud_s2_11_4,
)
from ._rabinowitz_richter import (
    rabinowitz_richter_1 as stroud_s2_9_2,
    rabinowitz_richter_2 as stroud_s2_9_4,
    rabinowitz_richter_4 as stroud_s2_11_3,
    rabinowitz_richter_5 as stroud_s2_13_1,
)
from ._radon import radon

from ._helpers import DiskScheme
from ..helpers import z, untangle, book


_citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_s2_5_1(symbolic=False):
    return radon(0, symbolic)


def stroud_s2_9_3(symbolic=False):
    # spherical product gauss 9
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pm_ = numpy.array([+1, -1])
    cos = numpy.vectorize(sympy.cos) if symbolic else numpy.cos
    sin = numpy.vectorize(sympy.sin) if symbolic else numpy.sin
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi

    r1, r2 = sqrt((6 - pm_ * sqrt(6)) / 10)

    a = (numpy.arange(10) + 1) * pi / 5
    x = numpy.array([cos(a), sin(a)]).T

    B0 = frac(1, 9)
    B1, B2 = (16 + pm_ * sqrt(6)) / 360

    data = [(B0, z(2)), (B1, r1 * x), (B2, r2 * x)]
    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Stroud S2 9-3", weights, points, 9, _citation)


__all__ = [
    "stroud_s2_3_1",
    "stroud_s2_3_2",
    "stroud_s2_4_1",
    "stroud_s2_5_1",
    "stroud_s2_5_2",
    "stroud_s2_7_1",
    "stroud_s2_7_2",
    "stroud_s2_9_1",
    "stroud_s2_9_2",
    "stroud_s2_9_3",
    "stroud_s2_9_4",
    "stroud_s2_9_5",
    "stroud_s2_11_1",
    "stroud_s2_11_2",
    "stroud_s2_11_3",
    "stroud_s2_11_4",
    "stroud_s2_13_1",
    "stroud_s2_13_2",
    "stroud_s2_15_1",
    "stroud_s2_15_2",
    "stroud_s2_17_1",
]
