# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .stenger import Stenger
from .stroud1967a import Stroud1967a
from .stroud1967b import Stroud1967b
from .stroud_secrest import StroudSecrest

from ..helpers import untangle, fsd, pm, pm_array0


class Stroud(object):
    """
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    """

    def __init__(self, n, index, symbolic=False):
        self.name = "Stroud_Enr2({})".format(index)
        self.dim = n

        scheme = {
            "3-1": lambda: StroudSecrest(n, "I", symbolic),
            "3-2": lambda: StroudSecrest(n, "III", symbolic),
            "5-1a": lambda: Stroud1967a(n, "a"),
            "5-1b": lambda: Stroud1967a(n, "b"),
            "5-2": lambda: StroudSecrest(n, "IV", symbolic),
            "5-3": lambda: Enr2_5_3(n, symbolic),
            "5-4": lambda: SphericalProductLobatto(n, symbolic),
            "5-5a": lambda: Enr2_5_5(n, True, symbolic),
            "5-5b": lambda: Enr2_5_5(n, False, symbolic),
            "5-6": lambda: Enr2_5_6(n, symbolic),
            "7-1a": lambda: Stroud1967b("2a", n, symbolic),
            "7-1b": lambda: Stroud1967b("2b", n, symbolic),
            "7-2": lambda: Stroud1967b("4", n, symbolic),
            "7-3a": lambda: Stenger(n, 7, "a"),
            "7-3b": lambda: Stenger(n, 7, "b"),
            "9-1a": lambda: Stenger(n, 9, "a"),
            "9-1b": lambda: Stenger(n, 9, "b"),
            "11-1a": lambda: Stenger(n, 11, "a"),
            "11-1b": lambda: Stenger(n, 11, "b"),
        }[index]()

        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return


class Enr2_5_3(object):
    def __init__(self, n, symbolic):
        sqrt = sympy.sqrt if symbolic else numpy.sqrt
        pi = sympy.pi if symbolic else numpy.pi
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.degree = 5

        r = sqrt(frac(n + 2, 4))
        s = sqrt(frac(n + 2, 2 * (n - 2)))
        A = frac(4, (n + 2) ** 2)
        B = frac((n - 2) ** 2, 2 ** n * (n + 2) ** 2)

        data = [(A, fsd(n, (r, 1))), (B, pm(n, s))]
        self.points, self.weights = untangle(data)
        self.weights *= sqrt(pi) ** n
        return


class SphericalProductLobatto(object):
    def __init__(self, n, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        pi = sympy.pi if symbolic else numpy.pi
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 5

        B0 = frac(2, (n + 2))
        data = [(B0, [n * [0]])]
        for k in range(1, n + 1):
            rk = sqrt(frac(k + 2, 2))
            s = sqrt(frac(1, 2))
            arr = [rk] + (n - k) * [s]
            idx = list(range(k - 1, n))
            alpha = frac(2 ** (k - n), (k + 1) * (k + 2))
            data += [(alpha, pm_array0(n, arr, idx))]

        self.points, self.weights = untangle(data)
        self.weights *= sqrt(pi) ** n

        return


class Enr2_5_5(object):
    def __init__(self, n, pos, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        pi = sympy.pi if symbolic else numpy.pi
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 5

        p_m = +1 if pos else -1

        # r is complex-valued for n >= 3
        r = sqrt((n + 2 + p_m * (n - 1) * sqrt(2 * (n + 2))) / (2 * n))
        s = sqrt((n + 2 - p_m * sqrt(2 * (n + 2))) / (2 * n))
        A = frac(2, n + 2)
        B = frac(1, 2 ** n * (n + 2))

        data = [(A, [n * [0]]), (B, fsd(n, (r, 1), (s, n - 1)))]

        self.points, self.weights = untangle(data)
        self.weights *= sqrt(pi) ** n

        return


class Enr2_5_6(object):
    def __init__(self, n, symbolic):
        assert n >= 5

        frac = sympy.Rational if symbolic else lambda x, y: x / y
        pi = sympy.pi if symbolic else numpy.pi
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 5

        sqrt2 = sqrt(2)
        sqrt2n1 = sqrt(2 * (n + 1))
        r = sqrt((n - sqrt2 + (n - 1) * sqrt2n1) / (2 * n))
        s = sqrt((n - sqrt2 - sqrt2n1) / (2 * n))
        t = sqrt((1 + sqrt2) / 2)
        A = frac(1, 2 ** n * (n + 1))

        data = [(A, fsd(n, (r, 1), (s, n - 1))), (A, pm(n, t))]

        self.points, self.weights = untangle(data)
        self.weights *= sqrt(pi) ** n
        return
