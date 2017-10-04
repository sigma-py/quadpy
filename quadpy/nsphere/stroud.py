# -*- coding: utf-8 -*-
#
from __future__ import division

from sympy import sqrt, Rational as fr

from ..helpers import untangle, fsd, pm_array0, pm

from .stroud1967 import Stroud1967
from .stroud1969 import Stroud1969
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
                (fr(1, 2*n), fsd(n, (1, 1))),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 3-2':
            self.degree = 3
            data = [
                (fr(1, 2**n), pm(n, sqrt(fr(1, n)))),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 5-1':
            self.degree = 5

            B1 = fr(4-n, 2*n*(n+2))
            B2 = fr(1, n * (n+2))

            data = [
                (B1, fsd(n, (1, 1))),
                (B2, fsd(n, (sqrt(fr(1, 2)), 2))),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 5-2':
            self.degree = 5

            B1 = fr(1, n * (n+2))
            B2 = fr(n, 2**n * (n+2))

            data = [
                (B1, fsd(n, (1, 1))),
                (B2, pm(n, sqrt(fr(1, n)))),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 5-3':
            self.degree = 5

            s = sqrt(fr(1, n+2))
            B = [
                fr(2**(k-n) * (n+2), n * (k+1) * (k+2))
                for k in range(1, n+1)
                ]
            r = [
                sqrt(fr(k+2, n+2))
                for k in range(1, n+1)
                ]
            data = [
                (B[k], pm_array0(n, [r[k]] + (n-k-1)*[s], range(k, n)))
                for k in range(n)
                ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 5-4':
            self.degree = 5

            s = sqrt(2*(n+2))
            u = sqrt((n + 2 + (n-1)*s) / n / (n+2))
            v = sqrt((n + 2 - s) / n / (n+2))

            data = [
                (fr(1, 2**n * n), fsd(n, (u, 1), (v, n-1))),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        elif index == 'Un 7-1':
            self.set_data(Stroud1967(n))
        elif index == 'Un 7-2':
            self.degree = 7

            A = fr(-n**2, 2**(n+3) * (n+2))
            B = fr((n+4)**2, 2**(n+3) * n * (n+2))

            r = sqrt(fr(1, n))
            s = sqrt(fr(5, n+4))
            t = sqrt(fr(1, n+4))

            data = [
                (A, pm(n, r)),
                (B, fsd(n, (s, 1), (t, n-1))),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        else:
            assert index == 'Un 11-1'
            self.set_data(Stroud1969(n))

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
