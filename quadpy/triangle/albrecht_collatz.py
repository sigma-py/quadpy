# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import untangle2


class AlbrechtCollatz(object):
    """
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://doi.org/10.1002/zamm.19580380102>
    """

    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.degree = 3

        self.data = {"s2": [[frac(2, 30), frac(1, 2)], [frac(9, 15), frac(1, 6)]]}
        self.bary, self.weights = untangle2(self.data)

        self.points = self.bary[:, 1:]
        self.weights /= 2
        return
