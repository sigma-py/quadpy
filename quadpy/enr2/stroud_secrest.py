# -*- coding: utf-8 -*-
#
from __future__ import division

from math import sqrt, pi

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
                (1/(n+1), sqrt(0.5) * _nsimplex(n))
                ]
        elif index == 'II':
            self.degree = 3
            nu = sqrt(n/2)
            data = [
                (1/(2*n), fsd(n, (nu, 1)))
                ]
        elif index == 'III':
            self.degree = 3
            nu = sqrt(0.5)
            data = [
                (1/(2**n), pm(n, nu))
                ]
        else:
            assert index == 'IV'
            self.degree = 5

            nu = sqrt((n+2) / 2)
            xi = sqrt((n+2) / 4)
            A = 2 / (n+2)
            B = (4-n) / 2 / (n+2)**2
            C = 1 / (n+2)**2

            data = [
                (A, numpy.array([numpy.full(n, 0.0)])),
                (B, fsd(n, (nu, 1))),
                (C, fsd(n, (xi, 2))),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= sqrt(pi)**n
        return


def _nsimplex(n):
    # construct the regular n-simplex points with 0 center
    return numpy.array([
        numpy.concatenate([
            -numpy.sqrt(
                (n + 1) / (n+1-numpy.arange(i)) / (n-numpy.arange(i))
                ),
            [numpy.sqrt((n+1) * (n-i) / (n+1-i))],
            numpy.zeros(n-i-1)
            ])
        for i in range(n)
        ] + [
        -numpy.sqrt(
            (n + 1) / (n+1-numpy.arange(n)) / (n-numpy.arange(n))
            )
        ])
