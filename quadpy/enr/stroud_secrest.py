# -*- coding: utf-8 -*-
#
"""
A.H. Stroud and D. Secrest,
Approximate integration formulas for certain spherically symmetric regions,
Math. Comp. 17 (1963), 105-135,
<https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
"""
from __future__ import division

import numpy
import scipy.special
import sympy

from ..helpers import untangle, fsd, pm
from ..enr2.stroud_secrest import _nsimplex


def i(n, symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    data = [(frac(1, n + 1), sqrt(n + 1) * _nsimplex(n, symbolic=symbolic))]
    return 2, data


def ii(n, symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    nu = sqrt(n * (n + 1))
    data = [(frac(1, 2 * n), fsd(n, (nu, 1)))]
    return 3, data


def iii(n, symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    nu = sqrt(n + 1)
    data = [(frac(1, 2 ** n), pm(n, nu))]
    return 3, data


def iv(n, symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    nu = sqrt((n + 2) * (n + 3))
    xi = sqrt(frac((n + 2) * (n + 3), 2))
    A = frac(2 * (2 * n + 3), (n + 2) * (n + 3))
    B = frac((4 - n) * (n + 1), 2 * (n + 2) ** 2 * (n + 3))
    C = frac(n + 1, (n + 2) ** 2 * (n + 3))

    data = [(A, numpy.full((1, n), 0)), (B, fsd(n, (nu, 1))), (C, fsd(n, (xi, 2)))]
    return 5, data


_gen = {"I": i, "II": ii, "III": iii, "IV": iv}


class StroudSecrest(object):
    keys = _gen.keys()

    def __init__(self, n, key, symbolic=False):
        self.dim = n
        self.degree, data = _gen[key](n, symbolic)
        self.points, self.weights = untangle(data)

        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt
        pi = sympy.pi if symbolic else numpy.pi
        gamma = sympy.gamma if symbolic else scipy.special.gamma

        self.weights *= 2 * sqrt(pi) ** n * gamma(n) / gamma(frac(n, 2))
        return
