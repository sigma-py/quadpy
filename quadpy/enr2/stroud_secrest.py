# -*- coding: utf-8 -*-
#
from __future__ import division

from sympy import sqrt, pi, Rational as fr

import numpy

from ..helpers import untangle, pm, fsd


class StroudSecrest(object):
    '''
    A.H. Stroud and D. Secrest,
    Approximate integration formulas for certain spherically symmetric regions,
    Math. Comp. 17 (1963), 105-135,
    <https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
    '''
    def __init__(self, n, index):
        self.dim = n
        if index == 'I':
            self.degree = 2
            data = [
                (fr(1, n+1), sqrt(fr(1, 2)) * _nsimplex(n))
                ]
        elif index == 'II':
            self.degree = 3
            nu = sqrt(fr(n, 2))
            data = [
                (fr(1, 2*n), fsd(n, (nu, 1)))
                ]
        elif index == 'III':
            self.degree = 3
            nu = sqrt(fr(1, 2))
            data = [
                (fr(1, 2**n), pm(n, nu))
                ]
        else:
            assert index == 'IV'
            self.degree = 5

            nu = sqrt(fr(n+2, 2))
            xi = sqrt(fr(n+2, 4))
            A = fr(2, n+2)
            B = fr(4-n, 2 * (n+2)**2)
            C = fr(1, (n+2)**2)

            data = [
                (A, numpy.full((1, n), 0)),
                (B, fsd(n, (nu, 1))),
                (C, fsd(n, (xi, 2))),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= sqrt(pi)**n
        return


def _nsimplex(n):
    # construct the regular n-simplex points with 0 center
    return numpy.array([
        [-sqrt(fr(n+1, (n+1-k) * (n-k))) for k in range(i)]
        + [sqrt(fr((n+1) * (n-i), n+1-i))]
        + (n-i-1) * [0]
        for i in range(n)
        ]
        + [
            [-sqrt(fr(n+1, (n+1-i) * (n-i))) for i in range(n)]
        ]
        )
