# -*- coding: utf-8 -*-
#
'''
A.H. Stroud and D. Secrest,
Approximate integration formulas for certain spherically symmetric regions,
Math. Comp. 17 (1963), 105-135,
<https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
'''
from __future__ import division

from sympy import sqrt, pi, gamma, Rational as fr

import numpy

from ..helpers import untangle, fsd, pm
from ..enr2.stroud_secrest import _nsimplex


def i(n):
    data = [
        (fr(1, n+1), sqrt(n+1) * _nsimplex(n, symbolic=True))
        ]
    return 2, data


def ii(n):
    nu = sqrt(n * (n+1))
    data = [
        (fr(1, 2*n), fsd(n, (nu, 1))),
        ]
    return 3, data


def iii(n):
    nu = sqrt(n+1)
    data = [
        (fr(1, 2**n), pm(n, nu)),
        ]
    return 3, data


def iv(n):
    nu = sqrt((n+2) * (n+3))
    xi = sqrt(fr((n+2) * (n+3), 2))
    A = fr(2*(2*n+3), (n+2) * (n+3))
    B = fr((4-n)*(n+1), 2 * (n+2)**2 * (n+3))
    C = fr(n+1, (n+2)**2 * (n+3))

    data = [
        (A, numpy.full((1, n), 0)),
        (B, fsd(n, (nu, 1))),
        (C, fsd(n, (xi, 2))),
        ]
    return 5, data


_gen = {
    'I': i,
    'II': ii,
    'III': iii,
    'IV': iv,
    }


class StroudSecrest(object):
    keys = _gen.keys()

    def __init__(self, n, key):
        self.dim = n
        self.degree, data = _gen[key](n)
        self.points, self.weights = untangle(data)
        self.weights *= 2 * sqrt(pi)**n * gamma(n) / gamma(fr(n, 2))
        return
