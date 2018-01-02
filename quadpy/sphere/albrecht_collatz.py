# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import cartesian_to_spherical_sympy
from ..helpers import untangle, pm_array0, fsd, pm


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self, index, symbolic=True):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        if index == 1:
            self.degree = 5
            r, s = [sqrt((5 + t * sqrt(5)) / 10) for t in [+1, -1]]
            data = [
                (frac(1, 12), pm_array0(3, [r, s], [0, 1])),
                (frac(1, 12), pm_array0(3, [r, s], [1, 2])),
                (frac(1, 12), pm_array0(3, [r, s], [2, 0])),
                ]
        elif index == 2:
            self.degree = 5
            r = 1
            s = sqrt(frac(1, 3))
            data = [
                (frac(8, 120), fsd(3, (r, 1))),
                (frac(9, 120), pm(3, s)),
                ]
        elif index == 3:
            self.degree = 5
            r = 1
            s = sqrt(frac(1, 2))
            data = [
                (frac(1, 30), fsd(3, (r, 1))),
                (frac(2, 30), fsd(3, (s, 2))),
                ]
        elif index == 4:
            self.degree = 5
            r, s = [sqrt((3 + t * sqrt(5)) / 6) for t in [+1, -1]]
            t = sqrt(frac(1, 3))
            data = [
                (frac(1, 20), pm_array0(3, [r, s], [0, 1])),
                (frac(1, 20), pm_array0(3, [r, s], [1, 2])),
                (frac(1, 20), pm_array0(3, [r, s], [2, 0])),
                (frac(1, 20), pm(3, t)),
                ]
        else:
            assert index == 5
            self.degree = 7

            r = 1
            s = sqrt(frac(1, 2))
            t = sqrt(frac(1, 3))

            data = [
                (frac(40, 840), fsd(3, (r, 1))),
                (frac(32, 840), fsd(3, (s, 2))),
                (frac(27, 840), pm(3, t)),
                ]

        self.points, self.weights = untangle(data)
        self.azimuthal_polar = cartesian_to_spherical_sympy(self.points)
        return
