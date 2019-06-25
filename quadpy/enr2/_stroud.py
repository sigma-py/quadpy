# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ._stenger import (
    stenger_7a as stroud_enr2_7_3a,
    stenger_7b as stroud_enr2_7_3b,
    stenger_9a as stroud_enr2_9_1a,
    stenger_9b as stroud_enr2_9_1b,
    stenger_11a as stroud_enr2_11_1a,
    stenger_11b as stroud_enr2_11_1b,
)
from ._stroud_1967_5 import (
    stroud_1967_5_a as stroud_enr2_5_1a,
    stroud_1967_5_b as stroud_enr2_5_1b,
)
from ._stroud_1967_7 import (
    stroud_1967_7_2a as stroud_enr2_7_1a,
    stroud_1967_7_2b as stroud_enr2_7_1b,
    stroud_1967_7_4 as stroud_enr2_7_2,
)
from ._stroud_secrest import (
    stroud_secrest_i as stroud_enr2_3_1,
    stroud_secrest_iii as stroud_enr2_3_2,
    stroud_secrest_iv as stroud_enr2_5_2,
)

from ._helpers import Enr2Scheme
from ..helpers import untangle, fsd, pm, pm_array0, book

citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_enr2_5_3(n, symbolic=False):
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    assert n > 2

    r = sqrt(frac(n + 2, 4))
    s = sqrt(frac(n + 2, 2 * (n - 2)))
    A = frac(4, (n + 2) ** 2)
    B = frac((n - 2) ** 2, 2 ** n * (n + 2) ** 2)

    data = [(A, fsd(n, (r, 1))), (B, pm(n, s))]
    points, weights = untangle(data)
    weights *= sqrt(pi) ** n
    return Enr2Scheme("Stroud Enr2 5-3", n, weights, points, 5, citation)


def stroud_enr2_5_4(n, symbolic=False):
    # Spherical product Lobatto
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    B0 = frac(2, (n + 2))
    data = [(B0, [n * [0]])]
    for k in range(1, n + 1):
        rk = sqrt(frac(k + 2, 2))
        s = sqrt(frac(1, 2))
        arr = [rk] + (n - k) * [s]
        idx = list(range(k - 1, n))
        alpha = frac(2 ** (k - n), (k + 1) * (k + 2))
        data += [(alpha, pm_array0(n, arr, idx))]

    points, weights = untangle(data)
    weights *= sqrt(pi) ** n
    return Enr2Scheme("Stroud Enr2 5-4", n, weights, points, 5, citation)


def _stroud_5_5(n, variant_a, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    p_m = +1 if variant_a else -1

    # r is complex-valued for n >= 3
    r = sqrt((n + 2 + p_m * (n - 1) * sqrt(2 * (n + 2))) / (2 * n))
    s = sqrt((n + 2 - p_m * sqrt(2 * (n + 2))) / (2 * n))
    A = frac(2, n + 2)
    B = frac(1, 2 ** n * (n + 2))

    data = [(A, [n * [0]]), (B, fsd(n, (r, 1), (s, n - 1)))]

    points, weights = untangle(data)
    weights *= sqrt(pi) ** n
    name = "Stroud Enr2 5-5{}".format("a" if variant_a else "b")
    return Enr2Scheme(name, n, weights, points, 5, citation)


def stroud_enr2_5_5a(n, symbolic=False):
    return _stroud_5_5(n, True, symbolic)


def stroud_enr2_5_5b(n, symbolic=False):
    return _stroud_5_5(n, False, symbolic)


def stroud_enr2_5_6(n, symbolic=False):
    assert n >= 5

    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    sqrt2 = sqrt(2)
    sqrt2n1 = sqrt(2 * (n + 1))
    r = sqrt((n - sqrt2 + (n - 1) * sqrt2n1) / (2 * n))
    s = sqrt((n - sqrt2 - sqrt2n1) / (2 * n))
    t = sqrt((1 + sqrt2) / 2)
    A = frac(1, 2 ** n * (n + 1))

    data = [(A, fsd(n, (r, 1), (s, n - 1))), (A, pm(n, t))]

    points, weights = untangle(data)
    weights *= sqrt(pi) ** n
    return Enr2Scheme("Stroud Enr2 5-6", n, weights, points, 5, citation)


__all__ = [
    "stroud_enr2_3_1",
    "stroud_enr2_3_2",
    "stroud_enr2_5_1a",
    "stroud_enr2_5_1b",
    "stroud_enr2_5_2",
    "stroud_enr2_5_3",
    "stroud_enr2_5_4",
    "stroud_enr2_5_5a",
    "stroud_enr2_5_5b",
    "stroud_enr2_5_6",
    "stroud_enr2_7_1a",
    "stroud_enr2_7_1b",
    "stroud_enr2_7_2",
    "stroud_enr2_7_3a",
    "stroud_enr2_7_3b",
    "stroud_enr2_9_1a",
    "stroud_enr2_9_1b",
    "stroud_enr2_11_1a",
    "stroud_enr2_11_1b",
]
