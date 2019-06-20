# -*- coding: utf-8 -*-
#
import numpy
import sympy

from ._helpers import cartesian_to_spherical, SphereScheme
from ..helpers import book

from ._albrecht_collatz import (
    albrecht_collatz_1 as stroud_u3_5_1,
    albrecht_collatz_2 as stroud_u3_5_2,
    albrecht_collatz_3 as stroud_u3_5_3,
    albrecht_collatz_4 as stroud_u3_5_4,
    albrecht_collatz_5 as stroud_u3_7_2,
)
from ._mclaren import (
    mclaren_01 as stroud_u3_3_1,
    mclaren_02 as stroud_u3_5_5,
    mclaren_03 as stroud_u3_7_1,
    mclaren_04 as stroud_u3_8_1,
    mclaren_05 as stroud_u3_9_1,
    mclaren_06 as stroud_u3_9_2,
    mclaren_07 as stroud_u3_9_3,
    mclaren_08 as stroud_u3_11_1,
    mclaren_09 as stroud_u3_11_3,
    mclaren_10 as stroud_u3_14_1,
)
from ..nsphere._stroud_1969 import stroud_1969

citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_u3_11_2(symbolic=False):
    scheme = stroud_1969(3, symbolic=symbolic)
    degree = scheme.degree
    weights = scheme.weights
    pi = sympy.pi if symbolic else numpy.pi
    weights /= 4 * pi
    points = scheme.points
    azimuthal_polar = cartesian_to_spherical(points)
    return SphereScheme(
        "Stroud U3 11-2", weights, points, azimuthal_polar, degree, citation
    )


__all__ = [
    "stroud_u3_3_1",
    "stroud_u3_5_1",
    "stroud_u3_5_2",
    "stroud_u3_5_3",
    "stroud_u3_5_4",
    "stroud_u3_5_5",
    "stroud_u3_7_2",
    "stroud_u3_7_1",
    "stroud_u3_8_1",
    "stroud_u3_9_1",
    "stroud_u3_9_2",
    "stroud_u3_9_3",
    "stroud_u3_11_1",
    "stroud_u3_11_2",
    "stroud_u3_11_3",
    "stroud_u3_14_1",
]
