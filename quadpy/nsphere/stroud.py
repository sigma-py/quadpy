# -*- coding: utf-8 -*-
#
from __future__ import division

import math

from ..helpers import untangle, fsd, pm_array0, pm

from .helpers import integrate_monomial_over_unit_nsphere


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    def __init__(self, n, index):
        self.dim = n
        if index == 'Un 3-1':
            self.degree = 3
            data = [
                (0.5/n, fsd(n, 1.0, 1)),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 3-2':
            self.degree = 3
            data = [
                (0.5**n, pm(n, math.sqrt(1.0/n))),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 5-1':
            self.degree = 5

            B1 = (4.0 - n) / (2.0*n*(n+2))
            B2 = 1.0 / n / (n+2)

            data = [
                (B1, fsd(n, 1.0, 1)),
                (B2, fsd(n, math.sqrt(0.5), 2)),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 5-2':
            self.degree = 5

            B1 = 1.0 / n / (n+2)
            B2 = n / 2**n / (n+2)

            data = [
                (B1, fsd(n, 1.0, 1)),
                (B2, pm(n, math.sqrt(1.0/n))),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 5-3':
            self.degree = 5

            s = math.sqrt(1.0 / (n+2))
            B = [
                2.0**(k-n) * (n+2) / n / (k+1) / (k+2)
                for k in range(1, n+1)
                ]
            r = [
                math.sqrt((k+2) / (n+2))
                for k in range(1, n+1)
                ]
            data = [
                (B[k], pm_array0(n, [r[k]] + (n-k-1)*[s], range(k, n)))
                for k in range(n)
                ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        else:
            assert False

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
