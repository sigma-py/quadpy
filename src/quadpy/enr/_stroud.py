import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import book, fsd, pm, untangle
from ._helpers import EnrScheme
from ._stroud_secrest import stroud_secrest_2 as stroud_enr_3_1
from ._stroud_secrest import stroud_secrest_3 as stroud_enr_3_2
from ._stroud_secrest import stroud_secrest_4 as stroud_enr_5_1

source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)

# ERR
# TODO find mistake
# def _gen5_2(n):
#     assert n != 3
#
#     r2 = -(n+1)*(n+3) + (n+3)*sqrt((n+1)*(2*n+3))
#     s2 = n*(n+1)*(n+3) - 2*(n+3)*sqrt((n+1)*(2*n+3)) / (n-3) / (n+2)
#     A = (n+1) * (n+4) / r2**2
#     B = (n+1) * (n+3) / 2**n / s2**2
#     r = sqrt(r2)
#     s = sqrt(s2)
#
#     data = [
#         (A, fsd(n, (r, 1))),
#         (B, pm(n, s)),
#         ]
#     return 5, data


def stroud_enr_5_3(n):
    """Spherical product Lobatto formula."""
    data = []
    s = sqrt(n + 3)
    for k in range(1, n + 1):
        rk = sqrt((k + 2) * (n + 3))
        Bk = frac(2 ** (k - n) * (n + 1), (k + 1) * (k + 2) * (n + 3))
        arr = (k - 1) * [0] + [rk] + (n - k) * [s]
        data += [(Bk, pm(arr))]
    B0 = 1 - sum(item[0] * len(item[1]) for item in data)
    data += [(B0, np.full((1, n), 0))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return EnrScheme("Stroud Enr 5-3", n, weights, points, 5, source)


def stroud_enr_5_4(n):
    r = sqrt(((n + 2) * (n + 3) + (n - 1) * (n + 3) * sqrt(2 * (n + 2))) / n)
    s = sqrt(((n + 2) * (n + 3) - (n + 3) * sqrt(2 * (n + 2))) / n)
    A = frac(4 * n + 6, (n + 2) * (n + 3))
    B = frac(n + 1, (n + 2) * (n + 3) * 2 ** n)
    data = [(A, np.full((1, n), 0)), (B, fsd(n, (r, 1), (s, n - 1)))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return EnrScheme("Stroud Enr 5-4", n, weights, points, 5, source, 1.684e-11)


# math domain error
# def _gen5_5(n):
#     r = sqrt((n*(n+1) - sqrt((n+1)*(4*n+6)) + (n-1)*(n+1)*sqrt(4*n+6)) / n)
#     s = sqrt((n*(n+1) - sqrt((n+1)*(4*n+6)) - (n+1)*sqrt(4*n+6)) / n)
#     t = n + 1 + sqrt((n+1)*(4*n+6))
#     A = 1 / (n+1) / 2**n
#     data = [
#         (A, fsd(n, (r, 1), (s, n-1))),
#         (A, pm(n, t)),
#         ]
#     return 5, data


# TODO find out what's wrong
# def _gen7_1(n):
#     assert 3 <= n <= 7
#
#     alpha = sqrt(3*(n+3)*(2*n+7)*(8-n))
#
#     # One could change some signs here and the formula would still work, but
#     # weights and points would be complex-valued.
#     r2 = (n+5) * (3*(n+3)*(8-n) - (n-2) * alpha) / (-2*n**2 + 6*n + 11)
#     s2 = (n+5) * (3*n*(n+3) - 2 * alpha) / (3*n**2 + 5*n - 56)
#     t2 = (n+5) * (6*(n+3) + 2 * alpha) / (2*n - 5)
#
#     B = (8-n) * (n+1) * (n+3) * (n+5) / r2**3
#     C = (n+1) * (n+3) * (n+5) / 2**n / s2**3
#     D = (n+1) * (n+3) * (n+5) / 2 / t2**3
#     A = 1.0 - 2*n*B - 2**n*C - 2*n*(n-1)*D
#
#     r = sqrt(r2)
#     s = sqrt(s2)
#     t = sqrt(t2)
#
#     data = [
#         (A, np.full((1, n), 0.0)),
#         (B, fsd(n, (r, 1))),
#         (C, pm(n, s)),
#         (D, fsd(n, (t, 2))),
#         ]
#     return 7, data


__all__ = [
    "stroud_enr_3_1",
    "stroud_enr_3_2",
    "stroud_enr_5_1",
    "stroud_enr_5_3",
    "stroud_enr_5_4",
]
