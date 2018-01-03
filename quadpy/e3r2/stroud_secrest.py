# -*- coding: utf-8 -*-
#
'''
A.H. Stroud and D. Secrest,
Approximate integration formulas for certain spherically symmetric regions,
Math. Comp. 17 (1963), 105-135,
<https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
'''
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, pm_roll, fsd, pm


def vii(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    # article:
    # nu, xi = numpy.sqrt((15 + plus_minus * 3*numpy.sqrt(5)))
    # A = 3/5
    # B = 1/30

    # book:
    nu, xi = [sqrt((5 - p_m * sqrt(5)) / 4) for p_m in [+1, -1]]
    A = frac(2, 5)
    B = frac(1, 20)

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, pm_roll(3, [nu, xi])),
        ]
    return 5, data


def viiia(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    r = sqrt(frac(5, 4))
    s = sqrt(frac(5, 2))
    data = [
        (frac(4, 25), fsd(3, (r, 1))),
        (frac(1, 200), pm(3, s)),
        ]
    return 5, data


def viiib(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    r = sqrt(frac(5, 2))
    s = sqrt(frac(5, 6))
    data = [
        (frac(2, 5), numpy.array([[0, 0, 0]])),
        (frac(1, 25), fsd(3, (r, 1))),
        (frac(9, 200), pm(3, s)),
        ]
    return 5, data


def ix(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    r, s = [sqrt((15 - p_m * 5*sqrt(5))/12) for p_m in [+1, -1]]
    t = sqrt(frac(5, 6))

    data = [
        (frac(2, 5), numpy.array([[0, 0, 0]])),
        (frac(3, 100), pm_roll(3, [r, s])),
        (frac(3, 100), pm(3, t)),
        ]
    return 5, data


def x(plus_minus, symbolic):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 7

    sqrt15 = sqrt(15)

    r = sqrt((15 + plus_minus * sqrt15) / 4)
    s = sqrt((6 - plus_minus * sqrt15) / 2)
    t = sqrt((9 + plus_minus * 2*sqrt15) / 2)
    A = (720 + plus_minus * 8*sqrt15) / 2205
    B = (270 - plus_minus * 46*sqrt15) / 15435
    C = (162 + plus_minus * 41*sqrt15) / 6174
    D = (783 - plus_minus * 202*sqrt15) / 24696

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, fsd(3, (r, 1))),
        (C, fsd(3, (s, 2))),
        (D, pm(3, t)),
        ]
    return degree, data


# pylint: disable=too-many-locals
def xi_(p_m, symbolic):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 7

    sqrt2 = sqrt(2)
    sqrt5 = sqrt(5)
    sqrt10 = sqrt(10)

    r = sqrt((25 + p_m * 15*sqrt2 + 5*sqrt5 + p_m * 3*sqrt10)/4)
    s = sqrt((25 + p_m * 15*sqrt2 - 5*sqrt5 - p_m * 3*sqrt10)/4)
    t = sqrt((3 - p_m * sqrt2) / 2)
    u = sqrt((9 - p_m * 3*sqrt2 - 3*sqrt5 + p_m * sqrt10) / 4)
    v = sqrt((9 - p_m * 3*sqrt2 + 3*sqrt5 - p_m * sqrt10) / 4)

    A = (80 + p_m * 8*sqrt2) / 245
    B = (395 - p_m * 279*sqrt2) / 13720
    C = (45 + p_m * 29*sqrt2) / 2744

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, pm_roll(3, [r, s])),
        (C, pm_roll(3, [u, v])),
        (C, pm(3, t)),
        ]
    return degree, data


_gen = {
    'VII': vii,
    'VIIIa': viiia,
    'VIIIb': viiib,
    'IX': ix,
    'Xa': lambda symbolic: x(+1, symbolic),
    'Xb': lambda symbolic: x(-1, symbolic),
    'XIa': lambda symbolic: xi_(+1, symbolic),
    'XIb': lambda symbolic: xi_(-1, symbolic),
    }


class StroudSecrest(object):
    keys = _gen.keys()

    def __init__(self, key, symbolic=False):
        self.name = 'StroudSecrest_E3r2({})'.format(key)

        self.degree, data = _gen[key](symbolic)
        self.points, self.weights = untangle(data)

        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        pi = sympy.pi if symbolic else numpy.pi

        self.weights *= sqrt(pi)**3
        return
