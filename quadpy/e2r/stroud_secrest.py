# -*- coding: utf-8 -*-
#
'''
A.H. Stroud and D. Secrest,
Approximate integration formulas for certain spherically symmetric regions,
Math. Comp. 17 (1963), 105-135,
<https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
'''
from __future__ import division

from math import sqrt, pi

import numpy

from ..helpers import untangle, pm_array, pm, fsd


def _v():
    nu = 2 * sqrt(5)
    xi = sqrt(5)
    eta = sqrt(15)

    data = [
        (0.7, numpy.array([[0.0, 0.0]])),
        (0.05, numpy.array([[+nu, 0], [-nu, 0]])),
        (0.05, pm_array([xi, eta])),
        ]
    return 5, data


def _vi():
    p_m = numpy.array([+1, -1])
    sqrt74255 = sqrt(74255)

    nu = sqrt(42)
    xi, eta = numpy.sqrt((6615 - p_m * 21 * sqrt74255) / 454)
    A = 5 / 588
    B, C = (5272105 + p_m * 18733 * sqrt74255) / 43661940

    data = [
        (A, fsd(2, (nu, 1))),
        (B, pm(2, xi)),
        (C, pm(2, eta)),
        ]
    return 7, data


_gen = {
    'V': _v,
    'VI': _vi,
    }


class StroudSecrest(object):
    keys = _gen.keys()

    def __init__(self, key):
        self.degree, data = _gen[key]()
        self.points, self.weights = untangle(data)
        self.weights *= 2 * pi
        return
