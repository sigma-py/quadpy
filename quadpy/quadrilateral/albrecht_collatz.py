# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _symm_r_0, _symm_s, _z, _pm, _pm2, unroll
from ..helpers import untangle


class AlbrechtCollatz(object):
    """
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://doi.org/10.1002/zamm.19580380102>
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "AlbrechtCollatz({})".format(index)
        if index == 1:
            self.degree = 3
            data = {
                "zero": [[frac(5, 12)]],
                "symm_r0": [[frac(1, 8), 1]],
                "symm_s": [[frac(1, 48), 1]],
            }
            self.points, self.weights = unroll(data)
        elif index == 2:
            self.degree = 5
            r = sqrt(frac(3, 5))
            s = sqrt(frac(1, 3))
            t = sqrt(frac(14, 15))
            data = {
                "zero": [[frac(2, 7)]],
                "pm": [[frac(5, 63), 0, t]],
                "pm2": [[frac(5, 36), r, s]],
            }
            self.points, self.weights = unroll(data)
        elif index == 3:
            self.degree = 5
            r = sqrt(frac(7, 15))
            s, t = [sqrt((7 + i * sqrt(24)) / 15) for i in [+1, -1]]
            data = {
                "zero": [[frac(2, 7)]],
                "pm": [
                    [frac(25, 168), r, r],
                    [frac(5, 48), +s, -t],
                    [frac(5, 48), +t, -s],
                ],
            }
            self.points, self.weights = unroll(data)
        else:
            assert index == 4
            self.degree = 5
            data = {
                "zero": [[frac(2, 45)]],
                "symm_r0": [[frac(2, 45), 1]],
                "symm_s": [[frac(1, 60), 1], [frac(8, 45), frac(1, 2)]],
            }
            self.points, self.weights = unroll(data)

        self.weights *= 4
        return
