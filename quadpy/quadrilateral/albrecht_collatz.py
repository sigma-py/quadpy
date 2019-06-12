# -*- coding: utf-8 -*-
#
"""
J. Albrecht, L. Collatz,
Zur numerischen Auswertung mehrdimensionaler Integrale,
ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
<https://doi.org/10.1002/zamm.19580380102>
"""
from __future__ import division

import numpy
import sympy

from .helpers import concat, symm_r0, symm_s, pmy, pm2, pm


class AlbrechtCollatz1(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.name = "AlbrechtCollatz(1)"
        self.degree = 3

        self.points, self.weights = concat(
            [
                ([[0, 0]], [frac(5, 12)]),
                symm_r0([[frac(1, 8), 1]]),
                symm_s([[frac(1, 48), 1]]),
            ]
        )
        self.weights *= 4
        return


class AlbrechtCollatz2(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "AlbrechtCollatz(2)"
        self.degree = 5
        r = sqrt(frac(3, 5))
        s = sqrt(frac(1, 3))
        t = sqrt(frac(14, 15))
        self.points, self.weights = concat(
            [
                ([[0, 0]], [frac(2, 7)]),
                pmy([[frac(5, 63), t]]),
                pm2([[frac(5, 36), r, s]]),
            ]
        )
        self.weights *= 4
        return


class AlbrechtCollatz3(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "AlbrechtCollatz(3)"
        self.degree = 5

        r = sqrt(frac(7, 15))
        s, t = [sqrt((7 + i * sqrt(24)) / 15) for i in [+1, -1]]
        self.points, self.weights = concat(
            [
                ([[0, 0]], [frac(2, 7)]),
                pm(
                    [
                        [frac(25, 168), r, r],
                        [frac(5, 48), +s, -t],
                        [frac(5, 48), +t, -s],
                    ]
                ),
            ]
        )
        self.weights *= 4
        return


class AlbrechtCollatz4(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.degree = 5

        self.points, self.weights = concat(
            [
                ([[0, 0]], [frac(2, 45)]),
                symm_r0([[frac(2, 45), 1]]),
                symm_s([[frac(1, 60), 1], [frac(8, 45), frac(1, 2)]]),
            ]
        )

        self.weights *= 4
        return


AlbrechtCollatz = {
    1: AlbrechtCollatz1,
    2: AlbrechtCollatz2,
    3: AlbrechtCollatz3,
    4: AlbrechtCollatz4,
}
