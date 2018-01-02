# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, pm, fsd


class StroudSecrest(object):
    '''
    A.H. Stroud and D. Secrest,
    Approximate integration formulas for certain spherically symmetric regions,
    Math. Comp. 17 (1963), 105-135,
    <https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
    '''
    def __init__(self, n, index, symbolic=True):
        sqrt = sympy.sqrt if symbolic else numpy.sqrt
        pi = sympy.pi if symbolic else numpy.pi
        frac = sympy.Rational if symbolic else lambda x, y: x/y

        self.dim = n
        if index == 'I':
            self.degree = 2
            data = [
                (
                    frac(1, n+1),
                    sqrt(frac(1, 2)) * _nsimplex(n, symbolic=symbolic)
                )
                ]
        elif index == 'II':
            self.degree = 3
            nu = sqrt(frac(n, 2))
            data = [
                (frac(1, 2*n), fsd(n, (nu, 1)))
                ]
        elif index == 'III':
            self.degree = 3
            nu = sqrt(frac(1, 2))
            data = [
                (frac(1, 2**n), pm(n, nu))
                ]
        else:
            assert index == 'IV'
            self.degree = 5

            nu = sqrt(frac(n+2, 2))
            xi = sqrt(frac(n+2, 4))
            A = frac(2, n+2)
            B = frac(4-n, 2 * (n+2)**2)
            C = frac(1, (n+2)**2)

            data = [
                (A, numpy.full((1, n), 0)),
                (B, fsd(n, (nu, 1))),
                (C, fsd(n, (xi, 2))),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= sqrt(pi)**n
        return


def _nsimplex(n, symbolic):
    # construct the regular n-simplex points with 0 center
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    frac = sympy.Rational if symbolic else lambda x, y: x/y

    return numpy.array([
        [-sqrt(frac(n+1, (n+1-k) * (n-k))) for k in range(i)]
        + [sqrt(frac((n+1) * (n-i), n+1-i))]
        + (n-i-1) * [0]
        for i in range(n)
        ]
        + [
            [-sqrt(frac(n+1, (n+1-i) * (n-i))) for i in range(n)]
        ]
        )
