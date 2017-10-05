# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
from sympy import sqrt, Rational as fr, pi

from ..helpers import untangle, pm_array, pm_array0, fsd, pm


class StroudSecrest(object):
    '''
    A.H. Stroud and D. Secrest,
    Approximate integration formulas for certain spherically symmetric regions,
    Math. Comp. 17 (1963), 105-135,
    <https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
    '''
    def __init__(self, index):
        if index == 'V':
            self.degree = 5

            nu = sqrt(2)
            xi = nu / 2
            eta = sqrt(6) / 2
            A = fr(1, 2)
            B = fr(1, 12)

            data = [
                (A, numpy.array([[0, 0]])),
                (B, pm_array0(2, [nu], [0])),
                (B, pm_array([xi, eta])),
                ]
        else:
            assert index == 'VI'
            self.degree = 7

            sqrt5 = sqrt(5)
            nu = sqrt(3)
            xi, eta = [sqrt((9 - p_m * 3*sqrt5) / 8) for p_m in [+1, -1]]
            A = fr(1, 36)
            B, C = [(5 + p_m * 2*sqrt5) / 45 for p_m in [+1, -1]]

            data = [
                (A, fsd(2, (nu, 1))),
                (B, pm(2, xi)),
                (C, pm(2, eta)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
