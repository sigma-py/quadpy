# -*- coding: utf-8 -*-
#
from __future__ import division

from math import sqrt, pi, lgamma, exp

import numpy

from ..helpers import untangle, fsd, pm
from ..enr2.stroud_secrest import _nsimplex


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
                (1/(n+1), sqrt(n+1) * _nsimplex(n))
                ]
        elif index == 'II':
            self.degree = 3
            nu = sqrt(n * (n+1))
            data = [
                (1/(2*n), fsd(n, (nu, 1))),
                ]
        elif index == 'III':
            self.degree = 3
            nu = sqrt(n+1)
            data = [
                (1/2**n, pm(n, nu)),
                ]
        else:
            assert index == 'IV'
            self.degree = 5

            nu = sqrt((n+2) * (n+3))
            xi = sqrt((n+2) * (n+3) / 2)
            A = 2*(2*n+3) / (n+2) / (n+3)
            B = (4-n)*(n+1) / 2 / (n+2)**2 / (n+3)
            C = (n+1) / (n+2)**2 / (n+3)

            data = [
                (A, numpy.full((1, n), 0.0)),
                (B, fsd(n, (nu, 1))),
                (C, fsd(n, (xi, 2))),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 2 * sqrt(pi)**n * exp(lgamma(n) - lgamma(n/2))
        return
