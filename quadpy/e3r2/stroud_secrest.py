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

from ..helpers import untangle, pm_roll, fsd, pm


def vii():
    plus_minus = numpy.array([+1, -1])

    # article:
    # nu, xi = numpy.sqrt((15 + plus_minus * 3*numpy.sqrt(5)))
    # A = 3/5
    # B = 1/30

    # book:
    nu, xi = numpy.sqrt((5 - plus_minus * numpy.sqrt(5)) / 4)
    A = 2/5
    B = 1/20

    data = [
        (A, numpy.array([[0.0, 0.0, 0.0]])),
        (B, pm_roll(3, [nu, xi])),
        ]
    return 5, data


def viiia():
    r = numpy.sqrt(5/4)
    s = numpy.sqrt(5/2)
    data = [
        (4/25, fsd(3, (r, 1))),
        (1/200, pm(3, s)),
        ]
    return 5, data


def viiib():
    r = numpy.sqrt(5/2)
    s = numpy.sqrt(5/6)
    data = [
        (2/5, numpy.array([[0.0, 0.0, 0.0]])),
        (1/25, fsd(3, (r, 1))),
        (9/200, pm(3, s)),
        ]
    return 5, data


def ix():
    plus_minus = numpy.array([+1, -1])

    r, s = numpy.sqrt((15 - plus_minus * 5*numpy.sqrt(5))/12)
    t = numpy.sqrt(5/6)

    data = [
        (2/5, numpy.array([[0.0, 0.0, 0.0]])),
        (3/100, pm_roll(3, [r, s])),
        (3/100, pm(3, t)),
        ]
    return 5, data


def x(plus_minus):
    degree = 7

    sqrt15 = numpy.sqrt(15)

    r = numpy.sqrt((15 + plus_minus * sqrt15) / 4)
    s = numpy.sqrt((6 - plus_minus * sqrt15) / 2)
    t = numpy.sqrt((9 + plus_minus * 2*sqrt15) / 2)
    A = (720 + plus_minus * 8*sqrt15) / 2205
    B = (270 - plus_minus * 46*sqrt15) / 15435
    C = (162 + plus_minus * 41*sqrt15) / 6174
    D = (783 - plus_minus * 202*sqrt15) / 24696

    data = [
        (A, numpy.array([[0.0, 0.0, 0.0]])),
        (B, fsd(3, (r, 1))),
        (C, fsd(3, (s, 2))),
        (D, pm(3, t)),
        ]
    return degree, data


def xi_(p_m):
    degree = 7

    sqrt2 = numpy.sqrt(2)
    sqrt5 = numpy.sqrt(5)
    sqrt10 = numpy.sqrt(10)

    r = numpy.sqrt((25 + p_m * 15*sqrt2 + 5*sqrt5 + p_m * 3*sqrt10)/4)
    s = numpy.sqrt((25 + p_m * 15*sqrt2 - 5*sqrt5 - p_m * 3*sqrt10)/4)
    t = numpy.sqrt((3 - p_m * sqrt2) / 2)
    u = numpy.sqrt((9 - p_m * 3*sqrt2 - 3*sqrt5 + p_m * sqrt10) / 4)
    v = numpy.sqrt((9 - p_m * 3*sqrt2 + 3*sqrt5 - p_m * sqrt10) / 4)

    A = (80 + p_m * 8*sqrt2) / 245
    B = (395 - p_m * 279*sqrt2) / 13720
    C = (45 + p_m * 29*sqrt2) / 2744

    data = [
        (A, numpy.array([[0.0, 0.0, 0.0]])),
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
    'Xa': lambda: x(+1),
    'Xb': lambda: x(-1),
    'XIa': lambda: xi_(+1),
    'XIb': lambda: xi_(-1),
    }


class StroudSecrest(object):
    keys = _gen.keys()

    def __init__(self, key):
        self.degree, data = _gen[key]()
        self.points, self.weights = untangle(data)
        self.weights *= numpy.sqrt(numpy.pi)**3
        return
