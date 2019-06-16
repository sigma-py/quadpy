# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from __future__ import division

import numpy
import sympy

from .hammer_stroud import hammer_stroud_11n, hammer_stroud_12n
from .stenger import stenger_7b, stenger_9a, stenger_9b, stenger_11a, stenger_11b
from .stroud_1957 import Stroud_1957
from .stroud_1966 import stroud_1966_a, stroud_1966_b, stroud_1966_c, stroud_1966_d
from .stroud_1967_5 import stroud_1967_5_a, stroud_1967_5_b
from .stroud_1967_7 import stroud_1967_7_a, stroud_1967_7_b, stroud_1967_7_c

from .helpers import volume_unit_ball, NBallScheme
from ..helpers import pm, untangle


def stroud_3_2(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = sqrt(frac(1, n + 2))
    data = [(frac(1, 2 ** n), pm(n, r))]
    points, weights = untangle(data)
    weights *= volume_unit_ball(n, symbolic=symbolic)
    return NBallScheme("Stroud 3-2", n, 3, weights, points)


Stroud = {
    "Sn 2-1": Stroud_1957,
    "Sn 3-1": lambda n, symbolic=False: hammer_stroud_11n(n, 0, symbolic),
    "Sn 3-2": stroud_3_2,
    "Sn 5-1a": stroud_1967_5_a,
    "Sn 5-1b": stroud_1967_5_b,
    "Sn 5-2": lambda n, symbolic=False: hammer_stroud_12n(n, 0, symbolic),
    "Sn 5-3": stroud_1966_a,
    "Sn 5-4": stroud_1966_b,
    "Sn 5-5": stroud_1966_c,
    "Sn 5-6": stroud_1966_d,
    "Sn 7-1a": stroud_1967_7_a,
    "Sn 7-1b": stroud_1967_7_b,
    "Sn 7-2": stroud_1967_7_c,
    # TODO fix
    # "Sn 7-3a": lambda: Stenger(n, 7, "a"),
    "Sn 7-3b": stenger_7b,
    "Sn 9-1a": stenger_9a,
    "Sn 9-1b": stenger_9b,
    "Sn 11-1a": stenger_11a,
    "Sn 11-1b": stenger_11b,
}
