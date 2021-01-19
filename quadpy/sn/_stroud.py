import math

import numpy as np
import sympy

from ..helpers import book, pm, untangle
from ._hammer_stroud import hammer_stroud_11n, hammer_stroud_12n
from ._helpers import SnScheme
from ._stenger import stenger_7a as stroud_sn_7_3a
from ._stenger import stenger_7b as stroud_sn_7_3b
from ._stenger import stenger_9a as stroud_sn_9_1a
from ._stenger import stenger_9b as stroud_sn_9_1b
from ._stenger import stenger_11a as stroud_sn_11_1a
from ._stenger import stenger_11b as stroud_sn_11_1b
from ._stroud_1957 import stroud_1957 as stroud_sn_2_1
from ._stroud_1966 import stroud_1966_a as stroud_sn_5_3
from ._stroud_1966 import stroud_1966_b as stroud_sn_5_4
from ._stroud_1966 import stroud_1966_c as stroud_sn_5_5
from ._stroud_1966 import stroud_1966_d as stroud_sn_5_6
from ._stroud_1967_5 import stroud_1967_5_a as stroud_sn_5_1a
from ._stroud_1967_5 import stroud_1967_5_b as stroud_sn_5_1b
from ._stroud_1967_7 import stroud_1967_7_a as stroud_sn_7_1a
from ._stroud_1967_7 import stroud_1967_7_b as stroud_sn_7_1b
from ._stroud_1967_7 import stroud_1967_7_c as stroud_sn_7_2

source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_sn_3_1(n):
    return hammer_stroud_11n(n, 0)


def stroud_sn_3_2(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt

    r = sqrt(frac(1, n + 2))
    data = [(frac(1, 2 ** n), pm(n * [r]))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Stroud Sn 3-2", n, weights, points, 3, source)


def stroud_sn_5_2(n):
    return hammer_stroud_12n(n, 0)


__all__ = [
    "stroud_sn_2_1",
    "stroud_sn_3_1",
    "stroud_sn_3_2",
    "stroud_sn_5_1a",
    "stroud_sn_5_1b",
    "stroud_sn_5_2",
    "stroud_sn_5_3",
    "stroud_sn_5_4",
    "stroud_sn_5_5",
    "stroud_sn_5_6",
    "stroud_sn_7_1a",
    "stroud_sn_7_1b",
    "stroud_sn_7_2",
    "stroud_sn_7_3a",
    "stroud_sn_7_3b",
    "stroud_sn_9_1a",
    "stroud_sn_9_1b",
    "stroud_sn_11_1a",
    "stroud_sn_11_1b",
]
