# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from __future__ import division

import numpy
import scipy.special
import sympy

from . import stroud_secrest

from ..helpers import untangle, pm_array0, fsd

# ERR
# TODO find mistake
# def _gen5_2(n, symbolic):
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


def _gen5_3(n, symbolic):
    """Spherical product Lobatto formula.
    """
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    data = []
    s = sqrt(n + 3)
    for k in range(1, n + 1):
        rk = sqrt((k + 2) * (n + 3))
        Bk = frac(2 ** (k - n) * (n + 1), (k + 1) * (k + 2) * (n + 3))
        arr = [rk] + (n - k) * [s]
        data += [(Bk, pm_array0(n, arr, range(k - 1, n)))]
    B0 = 1 - sum([item[0] * len(item[1]) for item in data])
    data += [(B0, numpy.full((1, n), 0))]
    return 5, data


def _gen5_4(n, symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    r = sqrt(((n + 2) * (n + 3) + (n - 1) * (n + 3) * sqrt(2 * (n + 2))) / n)
    s = sqrt(((n + 2) * (n + 3) - (n + 3) * sqrt(2 * (n + 2))) / n)
    A = frac(4 * n + 6, (n + 2) * (n + 3))
    B = frac(n + 1, (n + 2) * (n + 3) * 2 ** n)
    data = [(A, numpy.full((1, n), 0)), (B, fsd(n, (r, 1), (s, n - 1)))]
    return 5, data


# math domain error
# def _gen5_5(n, symbolic):
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
# def _gen7_1(n, symbolic):
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
#         (A, numpy.full((1, n), 0.0)),
#         (B, fsd(n, (r, 1))),
#         (C, pm(n, s)),
#         (D, fsd(n, (t, 2))),
#         ]
#     return 7, data


_gen = {
    "3-1": stroud_secrest.ii,
    "3-2": stroud_secrest.iii,
    "5-1": stroud_secrest.iv,
    # '5-2': _gen5_2,
    "5-3": _gen5_3,
    "5-4": _gen5_4,
    # '5-5': _gen5_5,
    # '7-1': _gen7_1,
}


class Stroud(object):
    keys = _gen.keys()

    def __init__(self, n, key, symbolic=False):
        self.name = "Stround_Enr({})".format(key)
        self.dim = n
        self.degree, data = _gen[key](n, symbolic=symbolic)
        self.points, self.weights = untangle(data)

        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt
        pi = sympy.pi if symbolic else numpy.pi
        gamma = sympy.gamma if symbolic else scipy.special.gamma

        self.weights *= 2 * sqrt(pi) ** n * gamma(n) / gamma(frac(n, 2))
        return
