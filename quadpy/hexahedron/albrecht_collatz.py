# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import fs_r00, fs_rr0, z
from ..helpers import untangle


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y

        self.degree = 3
        data = [
            (frac(1, 4), z()),
            (frac(1, 12), fs_r00(1)),
            (frac(1, 48), fs_rr0(1)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8
        return
