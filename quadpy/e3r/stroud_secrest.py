# -*- coding: utf-8 -*-
#
'''
A.H. Stroud and D. Secrest,
Approximate integration formulas for certain spherically symmetric regions,
Math. Comp. 17 (1963), 105-135,
<https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
'''
from __future__ import division

import warnings

import numpy
import sympy

from ..helpers import untangle, pm, fsd, pm_roll


def vii(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    nu, xi = [sqrt(15 - p_m * 3 * sqrt(5)) for p_m in [+1, -1]]
    A = frac(3, 5)
    B = frac(1, 30)

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, pm_roll(3, [xi, nu])),
        ]
    return 5, data


def viii(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    nu = sqrt(30)
    eta = sqrt(10)
    A = frac(3, 5)
    B = frac(2, 75)
    C = frac(3, 100)

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, fsd(3, (nu, 1))),
        (C, pm(3, eta)),
        ]
    return 5, data


def ix(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    eta = sqrt(10)
    xi, nu = [sqrt(15 - p_m * 5 * sqrt(5)) for p_m in [+1, -1]]
    A = frac(3, 5)
    B = frac(1, 50)

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, pm(3, eta)),
        (B, pm_roll(3, [xi, nu])),
        ]
    return 5, data


def x(symbolic):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt130 = sqrt(130)

    nu = sqrt((720 - 24*sqrt130) / 11)
    xi = sqrt(288 + 24*sqrt130)
    eta = sqrt((-216 + 24*sqrt130) / 7)
    A = (5175 - 13*sqrt130) / 8820
    B = (3870 + 283*sqrt130) / 493920
    C = (3204 - 281*sqrt130) / 197568
    # ERR in Stroud's book: 917568 vs. 197568
    D = (4239 + 373*sqrt130) / 197568

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, fsd(3, (nu, 1))),
        (C, fsd(3, (xi, 2))),
        (D, pm(3, eta)),
        ]

    # ERR
    # TODO find out what's wrong
    warnings.warn('Stroud-Secrest\'s scheme X for E_3^r has degree 3, not 7.')
    return 3, data


def xi_(symbolic):
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt5 = sqrt(5)
    sqrt39 = sqrt(39)
    sqrt195 = sqrt(195)

    nu, xi = [
        sqrt(-50 + p_m*10*sqrt5 + 10*sqrt39 - p_m*2*sqrt195)
        for p_m in [+1, -1]
        ]
    eta = sqrt(36 + 4*sqrt39)
    mu, lmbda = [
        sqrt(54 + p_m * 18*sqrt5 + 6*sqrt39 + p_m * 2*sqrt195)
        for p_m in [+1, -1]
        ]
    A = (1725 - 26*sqrt39) / 2940
    B = (1065 + 171*sqrt39) / 54880
    C = (297 - 47*sqrt39) / 32928

    data = [
        (A, numpy.array([[0, 0, 0]])),
        (B, pm_roll(3, [xi, nu])),
        (C, pm(3, eta)),
        (C, pm_roll(3, [lmbda, mu])),
        ]
    return 7, data


_gen = {
    'VII': vii,
    'VIII': viii,
    'IX': ix,
    'X': x,
    'XI': xi_,
    }


class StroudSecrest(object):
    keys = _gen.keys()

    def __init__(self, key, symbolic=False):
        self.degree, data = _gen[key](symbolic)
        self.points, self.weights = untangle(data)
        pi = sympy.pi if symbolic else numpy.pi
        self.weights *= 8 * pi
        return
