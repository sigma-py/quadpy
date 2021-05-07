import math

import numpy as np
import sympy

from ..helpers import book, fsd, pm, untangle
from ._helpers import Enr2Scheme
from ._stenger import stenger_7a as stroud_enr2_7_3a
from ._stenger import stenger_7b as stroud_enr2_7_3b
from ._stenger import stenger_9a as stroud_enr2_9_1a
from ._stenger import stenger_9b as stroud_enr2_9_1b
from ._stenger import stenger_11a as stroud_enr2_11_1a
from ._stenger import stenger_11b as stroud_enr2_11_1b
from ._stroud_1967_5 import stroud_1967_5_a as stroud_enr2_5_1a
from ._stroud_1967_5 import stroud_1967_5_b as stroud_enr2_5_1b
from ._stroud_1967_7 import stroud_1967_7_2a as stroud_enr2_7_1a
from ._stroud_1967_7 import stroud_1967_7_2b as stroud_enr2_7_1b
from ._stroud_1967_7 import stroud_1967_7_4 as stroud_enr2_7_2
from ._stroud_secrest import stroud_secrest_1 as stroud_enr2_3_1
from ._stroud_secrest import stroud_secrest_2 as stroud_enr2_3_2
from ._stroud_secrest import stroud_secrest_4 as stroud_enr2_5_2

source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_enr2_5_3(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    assert n > 2

    r = sqrt(frac(n + 2, 4))
    s = sqrt(frac(n + 2, 2 * (n - 2)))
    A = frac(4, (n + 2) ** 2)
    B = frac((n - 2) ** 2, 2 ** n * (n + 2) ** 2)

    data = [(A, fsd(n, (r, 1))), (B, pm(n * [s]))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud Enr2 5-3", n, weights, points, 5, source)


def stroud_enr2_5_4(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    # Spherical product Lobatto
    B0 = frac(2, (n + 2))
    data = [(B0, [n * [0]])]
    for k in range(1, n + 1):
        rk = sqrt(frac(k + 2, 2))
        s = sqrt(frac(1, 2))
        arr = (k - 1) * [0] + [rk] + (n - k) * [s]
        alpha = frac(2 ** (k - n), (k + 1) * (k + 2))
        data += [(alpha, pm(arr))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud Enr2 5-4", n, weights, points, 5, source)


def _stroud_5_5(n, variant_a, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    p_m = +1 if variant_a else -1

    # r is complex-valued for n >= 3
    r = sqrt((n + 2 + p_m * (n - 1) * sqrt(2 * (n + 2))) / (2 * n))
    s = sqrt((n + 2 - p_m * sqrt(2 * (n + 2))) / (2 * n))
    A = frac(2, n + 2)
    B = frac(1, 2 ** n * (n + 2))

    data = [(A, [n * [0]]), (B, fsd(n, (r, 1), (s, n - 1)))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    variant = "a" if variant_a else "b"
    return Enr2Scheme(f"Stroud Enr2 5-5{variant}", n, weights, points, 5, source)


def stroud_enr2_5_5a(n):
    return _stroud_5_5(n, True)


def stroud_enr2_5_5b(n):
    return _stroud_5_5(n, False)


def stroud_enr2_5_6(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    assert n >= 5

    sqrt2 = sqrt(2)
    sqrt2n1 = sqrt(2 * (n + 1))
    r = sqrt((n - sqrt2 + (n - 1) * sqrt2n1) / (2 * n))
    s = sqrt((n - sqrt2 - sqrt2n1) / (2 * n))
    t = sqrt((1 + sqrt2) / 2)
    A = frac(1, 2 ** n * (n + 1))

    data = [(A, fsd(n, (r, 1), (s, n - 1))), (A, pm(n * [t]))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud Enr2 5-6", n, weights, points, 5, source)


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
