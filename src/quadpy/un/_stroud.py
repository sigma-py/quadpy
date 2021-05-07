from sympy import Rational as frac
from sympy import sqrt

from ..helpers import book, fsd, pm, untangle
from ._helpers import UnScheme
from ._stroud_1967 import stroud_1967 as stroud_un_7_1
from ._stroud_1969 import stroud_1969 as stroud_un_11_1

source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_un_3_1(n):
    degree = 3
    data = [(frac(1, 2 * n), fsd(n, (1, 1)))]
    points, weights = untangle(data)
    return UnScheme("Stroud Un 3-1", n, weights, points, degree, source)


def stroud_un_3_2(n):
    degree = 3
    data = [(frac(1, 2 ** n), pm(n * [sqrt(frac(1, n))]))]
    points, weights = untangle(data)
    return UnScheme("Stroud Un 3-2", n, weights, points, degree, source)


def stroud_un_5_1(n):
    degree = 5

    B1 = frac(4 - n, 2 * n * (n + 2))
    B2 = frac(1, n * (n + 2))

    data = [(B1, fsd(n, (1, 1))), (B2, fsd(n, (sqrt(frac(1, 2)), 2)))]

    points, weights = untangle(data)
    return UnScheme("Stroud Un 5-1", n, weights, points, degree, source)


def stroud_un_5_2(n):
    degree = 5

    B1 = frac(1, n * (n + 2))
    B2 = frac(n, 2 ** n * (n + 2))

    data = [(B1, fsd(n, (1, 1))), (B2, pm(n * [sqrt(frac(1, n))]))]

    points, weights = untangle(data)
    return UnScheme("Stroud Un 5-2", n, weights, points, degree, source)


def stroud_un_5_3(n):
    degree = 5

    s = sqrt(frac(1, n + 2))
    B = [frac(2 ** (k - n) * (n + 2), n * (k + 1) * (k + 2)) for k in range(1, n + 1)]
    r = [sqrt(frac(k + 2, n + 2)) for k in range(1, n + 1)]
    data = [(B[k], pm(k * [0] + [r[k]] + (n - k - 1) * [s])) for k in range(n)]

    points, weights = untangle(data)
    return UnScheme("Stroud Un 5-3", n, weights, points, degree, source)


def stroud_un_5_4(n):
    degree = 5

    s = sqrt(2 * (n + 2))
    u = sqrt((n + 2 + (n - 1) * s) / n / (n + 2))
    v = sqrt((n + 2 - s) / n / (n + 2))

    data = [(frac(1, 2 ** n * n), fsd(n, (u, 1), (v, n - 1)))]

    points, weights = untangle(data)
    return UnScheme("Stroud Un 5-4", n, weights, points, degree, source)


def stroud_un_7_2(n):
    degree = 7

    A = frac(-(n ** 2), 2 ** (n + 3) * (n + 2))
    B = frac((n + 4) ** 2, 2 ** (n + 3) * n * (n + 2))

    r = sqrt(frac(1, n))
    s = sqrt(frac(5, n + 4))
    t = sqrt(frac(1, n + 4))

    data = [(A, pm(n * [r])), (B, fsd(n, (s, 1), (t, n - 1)))]

    points, weights = untangle(data)
    return UnScheme("Stroud Un 7-2", n, weights, points, degree, source)


__all__ = [
    "stroud_un_3_1",
    "stroud_un_3_2",
    "stroud_un_5_1",
    "stroud_un_5_2",
    "stroud_un_5_3",
    "stroud_un_5_4",
    "stroud_un_7_1",
    "stroud_un_7_2",
    "stroud_un_11_1",
]
