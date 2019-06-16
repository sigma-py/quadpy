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

from .stenger import (
    stenger_7a,
    stenger_7b,
    stenger_9a,
    stenger_9b,
    stenger_11a,
    stenger_11b,
)
from .stroud_1967_5 import stroud_1967_5_a, stroud_1967_5_b
from .stroud_1967_7 import stroud_1967_7_2a, stroud_1967_7_2b, stroud_1967_7_4
from .stroud_secrest import stroud_secrest_i, stroud_secrest_iii, stroud_secrest_iv

from .helpers import Enr2Scheme
from ..helpers import untangle, fsd, pm, pm_array0


def stroud_5_3(n, symbolic=False):
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    r = sqrt(frac(n + 2, 4))
    s = sqrt(frac(n + 2, 2 * (n - 2)))
    A = frac(4, (n + 2) ** 2)
    B = frac((n - 2) ** 2, 2 ** n * (n + 2) ** 2)

    data = [(A, fsd(n, (r, 1))), (B, pm(n, s))]
    points, weights = untangle(data)
    weights *= sqrt(pi) ** n
    return Enr2Scheme("Stroud 5-3", n, 5, weights, points)


def stroud_5_4(n, symbolic=False):
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
    return Enr2Scheme("Stroud 5-4", n, 5, weights, points)


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
    name = "Stroud 5-5{}".format("a" if variant_a else "b")
    return Enr2Scheme(name, n, 5, weights, points)


def stroud_5_5a(n, symbolic=False):
    return _stroud_5_5(n, True, symbolic)


def stroud_5_5b(n, symbolic=False):
    return _stroud_5_5(n, False, symbolic)


def stroud_5_6(n, symbolic=False):
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
    return Enr2Scheme("Stroud 5-6", n, 5, weights, points)


Stroud = {
    "3-1": stroud_secrest_i,
    "3-2": stroud_secrest_iii,
    "5-1a": stroud_1967_5_a,
    "5-1b": stroud_1967_5_b,
    "5-2": stroud_secrest_iv,
    "5-3": stroud_5_3,
    "5-4": stroud_5_4,
    "5-5a": stroud_5_5a,
    "5-5b": stroud_5_5b,
    "5-6": stroud_5_6,
    "7-1a": stroud_1967_7_2a,
    "7-1b": stroud_1967_7_2b,
    "7-2": stroud_1967_7_4,
    "7-3a": stenger_7a,
    "7-3b": stenger_7b,
    "9-1a": stenger_9a,
    "9-1b": stenger_9b,
    "11-1a": stenger_11a,
    "11-1b": stenger_11b,
}
