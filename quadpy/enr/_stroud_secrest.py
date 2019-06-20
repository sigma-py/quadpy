# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import scipy.special
import sympy

from ._helpers import EnrScheme
from ..helpers import untangle, fsd, pm, article
from ..enr2.stroud_secrest import _nsimplex


citation = article(
    authors=["A.H. Stroud", "D. Secrest"],
    title="Approximate integration formulas for certain spherically symmetric regions",
    journal="Math. Comp.",
    volume="17",
    year="1963",
    pages="105-135",
    url="https://doi.org/10.1090/S0025-5718-1963-0161473-0",
)


def stroud_secrest_i(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    gamma = sympy.gamma if symbolic else scipy.special.gamma

    data = [(frac(1, n + 1), sqrt(n + 1) * _nsimplex(n, symbolic=symbolic))]
    points, weights = untangle(data)
    weights *= 2 * sqrt(pi) ** n * gamma(n) / gamma(frac(n, 2))
    return EnrScheme("Stroud-Secrest I", n, weights, points, 2, citation)


def stroud_secrest_ii(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    gamma = sympy.gamma if symbolic else scipy.special.gamma

    nu = sqrt(n * (n + 1))
    data = [(frac(1, 2 * n), fsd(n, (nu, 1)))]
    points, weights = untangle(data)
    weights *= 2 * sqrt(pi) ** n * gamma(n) / gamma(frac(n, 2))
    return EnrScheme("Stroud-Secrest II", n, weights, points, 3, citation)


def stroud_secrest_iii(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    gamma = sympy.gamma if symbolic else scipy.special.gamma

    nu = sqrt(n + 1)
    data = [(frac(1, 2 ** n), pm(n, nu))]
    points, weights = untangle(data)
    weights *= 2 * sqrt(pi) ** n * gamma(n) / gamma(frac(n, 2))
    return EnrScheme("Stroud-Secrest III", n, weights, points, 3, citation)


def stroud_secrest_iv(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    gamma = sympy.gamma if symbolic else scipy.special.gamma

    nu = sqrt((n + 2) * (n + 3))
    xi = sqrt(frac((n + 2) * (n + 3), 2))
    A = frac(2 * (2 * n + 3), (n + 2) * (n + 3))
    B = frac((4 - n) * (n + 1), 2 * (n + 2) ** 2 * (n + 3))
    C = frac(n + 1, (n + 2) ** 2 * (n + 3))

    data = [(A, numpy.full((1, n), 0)), (B, fsd(n, (nu, 1))), (C, fsd(n, (xi, 2)))]
    points, weights = untangle(data)
    weights *= 2 * sqrt(pi) ** n * gamma(n) / gamma(frac(n, 2))
    return EnrScheme("Stroud-Secrest IV", n, weights, points, 5, citation)
