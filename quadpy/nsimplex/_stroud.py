# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._hammer_stroud import (
    hammer_stroud_1a as stroud_tn_2_1a,
    hammer_stroud_1b as stroud_tn_2_1b,
    hammer_stroud_2 as stroud_tn_3_1,
)
from ._lauffer import (
    lauffer_1 as stroud_tn_1_2,
    lauffer_2 as stroud_tn_2_2,
    lauffer_3 as stroud_tn_3_9,
    lauffer_4 as stroud_tn_4_1,
    lauffer_5 as stroud_tn_5_2,
)
from ._stroud_1961 import stroud_1961 as stroud_tn_3_3
from ._stroud_1964 import stroud_1964a as stroud_tn_3_6a, stroud_1964b as stroud_tn_3_6b
from ._stroud_1966 import (
    stroud_1966_i as stroud_tn_3_2,
    stroud_1966_ii as stroud_tn_3_4,
    stroud_1966_iii as stroud_tn_3_5,
    stroud_1966_iv as stroud_tn_3_7,
    stroud_1966_v as stroud_tn_3_8,
    stroud_1966_vi as stroud_tn_3_10,
    stroud_1966_vii as stroud_tn_3_11,
)
from ._stroud_1969 import stroud_1969 as stroud_tn_5_1

from ..helpers import untangle, book
from ._helpers import NSimplexScheme


citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_tn_1_1(n, symbolic):
    # midpoint rule
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 1
    data = [(1, numpy.full((1, n + 1), frac(1, n + 1)))]
    bary, weights = untangle(data)
    return NSimplexScheme("Stroud Tn 1-1", n, weights, bary, degree, citation)


__all__ = [
    "stroud_tn_1_1",
    "stroud_tn_1_2",
    "stroud_tn_2_1a",
    "stroud_tn_2_1b",
    "stroud_tn_2_2",
    "stroud_tn_3_1",
    "stroud_tn_3_2",
    "stroud_tn_3_3",
    "stroud_tn_3_4",
    "stroud_tn_3_5",
    "stroud_tn_3_6a",
    "stroud_tn_3_6b",
    "stroud_tn_3_7",
    "stroud_tn_3_8",
    "stroud_tn_3_9",
    "stroud_tn_3_10",
    "stroud_tn_3_11",
    "stroud_tn_4_1",
    "stroud_tn_5_1",
    "stroud_tn_5_2",
]
