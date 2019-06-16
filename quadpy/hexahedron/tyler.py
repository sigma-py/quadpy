# -*- coding: utf-8 -*-
#
"""
G.W. Tyler,
Numerical integration of functions of several variables,
Canad. J. Math. 5(1953), 393-412,
<https://doi.org/10.4153/CJM-1953-044-1>.
"""
from __future__ import division

import sympy

from .helpers import fs_r00, pm_rrr, z, HexahedronScheme
from ..helpers import untangle


def tyler_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [(frac(1, 6), fs_r00(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Tyler 1", 3, weights, points)


def tyler_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [
        (-frac(62, 45), z()),
        (frac(16, 45), fs_r00(frac(1, 2))),
        (frac(1, 45), fs_r00(1)),
        (frac(1, 72), pm_rrr(1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Tyler 2", 5, weights, points)


Tyler = {1: tyler_1, 2: tyler_2}
