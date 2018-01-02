# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, pm


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        pi = sympy.pi if symbolic else numpy.pi

        self.name = 'AlbrechtCollatz'
        self.degree = 3
        data = [
            # ERR Wrongly stated in Stroud as sqrt(1/2) instead of 1/2
            (frac(1, 4), pm(2, frac(1, 2))),
            ]
        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
