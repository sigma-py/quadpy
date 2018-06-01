# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, fsd, pm_array0, pm

from .stroud1967 import Stroud1967
from .stroud1969 import Stroud1969
from .helpers import integrate_monomial_over_unit_nsphere


class Stroud(object):
    """
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    """

    def __init__(self, n, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "Stround Un({})".format(index)
        self.dim = n
        if index == "Un 3-1":
            self.degree = 3
            data = [(frac(1, 2 * n), fsd(n, (1, 1)))]
            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0], symbolic)
        elif index == "Un 3-2":
            self.degree = 3
            data = [(frac(1, 2 ** n), pm(n, sqrt(frac(1, n))))]
            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0], symbolic)
        elif index == "Un 5-1":
            self.degree = 5

            B1 = frac(4 - n, 2 * n * (n + 2))
            B2 = frac(1, n * (n + 2))

            data = [(B1, fsd(n, (1, 1))), (B2, fsd(n, (sqrt(frac(1, 2)), 2)))]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0], symbolic)
        elif index == "Un 5-2":
            self.degree = 5

            B1 = frac(1, n * (n + 2))
            B2 = frac(n, 2 ** n * (n + 2))

            data = [(B1, fsd(n, (1, 1))), (B2, pm(n, sqrt(frac(1, n))))]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0], symbolic)
        elif index == "Un 5-3":
            self.degree = 5

            s = sqrt(frac(1, n + 2))
            B = [
                frac(2 ** (k - n) * (n + 2), n * (k + 1) * (k + 2))
                for k in range(1, n + 1)
            ]
            r = [sqrt(frac(k + 2, n + 2)) for k in range(1, n + 1)]
            data = [
                (B[k], pm_array0(n, [r[k]] + (n - k - 1) * [s], range(k, n)))
                for k in range(n)
            ]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0], symbolic)
        elif index == "Un 5-4":
            self.degree = 5

            s = sqrt(2 * (n + 2))
            u = sqrt((n + 2 + (n - 1) * s) / n / (n + 2))
            v = sqrt((n + 2 - s) / n / (n + 2))

            data = [(frac(1, 2 ** n * n), fsd(n, (u, 1), (v, n - 1)))]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0], symbolic)
        elif index == "Un 7-1":
            self.set_data(Stroud1967(n, symbolic=symbolic))
        elif index == "Un 7-2":
            self.degree = 7

            A = frac(-n ** 2, 2 ** (n + 3) * (n + 2))
            B = frac((n + 4) ** 2, 2 ** (n + 3) * n * (n + 2))

            r = sqrt(frac(1, n))
            s = sqrt(frac(5, n + 4))
            t = sqrt(frac(1, n + 4))

            data = [(A, pm(n, r)), (B, fsd(n, (s, 1), (t, n - 1)))]

            self.points, self.weights = untangle(data)
            self.weights *= integrate_monomial_over_unit_nsphere(n * [0], symbolic)
        else:
            assert index == "Un 11-1"
            self.set_data(Stroud1969(n, symbolic=symbolic))

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
