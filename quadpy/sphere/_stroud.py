# -*- coding: utf-8 -*-
#
import sympy

from ..helpers import book
from ..nsphere._stroud_1969 import stroud_1969
from ._albrecht_collatz import albrecht_collatz_1 as stroud_u3_5_1
from ._albrecht_collatz import albrecht_collatz_2 as stroud_u3_5_2
from ._albrecht_collatz import albrecht_collatz_3 as stroud_u3_5_3
from ._albrecht_collatz import albrecht_collatz_4 as stroud_u3_5_4
from ._albrecht_collatz import albrecht_collatz_5 as stroud_u3_7_2
from ._helpers import SphereScheme, cartesian_to_spherical
from ._mclaren import mclaren_01 as stroud_u3_3_1
from ._mclaren import mclaren_02 as stroud_u3_5_5
from ._mclaren import mclaren_03 as stroud_u3_7_1
from ._mclaren import mclaren_04 as stroud_u3_8_1
from ._mclaren import mclaren_05 as stroud_u3_9_1
from ._mclaren import mclaren_06 as stroud_u3_9_2
from ._mclaren import mclaren_07 as stroud_u3_9_3
from ._mclaren import mclaren_08 as stroud_u3_11_1
from ._mclaren import mclaren_09 as stroud_u3_11_3
from ._mclaren import mclaren_10 as stroud_u3_14_1

citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_u3_11_2():
    scheme = stroud_1969(3)
    degree = scheme.degree
    weights = scheme.weights
    weights /= 4 * sympy.pi
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
