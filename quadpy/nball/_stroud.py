# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ._hammer_stroud import hammer_stroud_11n, hammer_stroud_12n
from ._stenger import (
    # TODO fix
    # stenger_7a as stroud_sn_7_3a,
    stenger_7b as stroud_sn_7_3b,
    stenger_9a as stroud_sn_9_1a,
    stenger_9b as stroud_sn_9_1b,
    stenger_11a as stroud_sn_11_1a,
    stenger_11b as stroud_sn_11_1b,
)
from ._stroud_1957 import stroud_1957 as stroud_sn_2_1
from ._stroud_1966 import (
    stroud_1966_a as stroud_sn_5_3,
    stroud_1966_b as stroud_sn_5_4,
    stroud_1966_c as stroud_sn_5_5,
    stroud_1966_d as stroud_sn_5_6,
)
from ._stroud_1967_5 import (
    stroud_1967_5_a as stroud_sn_5_1a,
    stroud_1967_5_b as stroud_sn_5_1b,
)
from ._stroud_1967_7 import (
    stroud_1967_7_a as stroud_sn_7_1a,
    stroud_1967_7_b as stroud_sn_7_1b,
    stroud_1967_7_c as stroud_sn_7_2,
)

from ._helpers import volume_unit_ball, NBallScheme
from ..helpers import pm, untangle, book

citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_sn_3_1(n, symbolic=False):
    return hammer_stroud_11n(n, 0, symbolic)


def stroud_sn_3_2(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = sqrt(frac(1, n + 2))
    data = [(frac(1, 2 ** n), pm(n, r))]
    points, weights = untangle(data)
    weights *= volume_unit_ball(n, symbolic=symbolic)
    return NBallScheme("Stroud Sn 3-2", n, weights, points, 3, citation)


def stroud_sn_5_2(n, symbolic=False):
    return hammer_stroud_12n(n, 0, symbolic)


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
    # "stroud_sn_7_3a",
    "stroud_sn_7_3b",
    "stroud_sn_9_1a",
    "stroud_sn_9_1b",
    "stroud_sn_11_1a",
    "stroud_sn_11_1b",
]
