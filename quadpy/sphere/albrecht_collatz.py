# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from .helpers import cartesian_to_spherical_sympy
from ..helpers import untangle, pm_array0, fsd, pm


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self, index):
        if index == 1:
            self.degree = 5
            r, s = [sqrt((5 + t * sqrt(5)) / 10) for t in [+1, -1]]
            data = [
                (fr(1, 12), pm_array0(3, [r, s], [0, 1])),
                (fr(1, 12), pm_array0(3, [r, s], [1, 2])),
                (fr(1, 12), pm_array0(3, [r, s], [2, 0])),
                ]
        elif index == 2:
            self.degree = 5
            r = 1
            s = sqrt(fr(1, 3))
            data = [
                (fr(8, 120), fsd(3, (r, 1))),
                (fr(9, 120), pm(3, s)),
                ]
        elif index == 3:
            self.degree = 5
            r = 1
            s = sqrt(fr(1, 2))
            data = [
                (fr(1, 30), fsd(3, (r, 1))),
                (fr(2, 30), fsd(3, (s, 2))),
                ]
        elif index == 4:
            self.degree = 5
            r, s = [sqrt((3 + t * sqrt(5)) / 6) for t in [+1, -1]]
            t = sqrt(fr(1, 3))
            data = [
                (fr(1, 20), pm_array0(3, [r, s], [0, 1])),
                (fr(1, 20), pm_array0(3, [r, s], [1, 2])),
                (fr(1, 20), pm_array0(3, [r, s], [2, 0])),
                (fr(1, 20), pm(3, t)),
                ]
        else:
            assert index == 5
            self.degree = 7

            r = 1
            s = sqrt(fr(1, 2))
            t = sqrt(fr(1, 3))

            data = [
                (fr(40, 840), fsd(3, (r, 1))),
                (fr(32, 840), fsd(3, (s, 2))),
                (fr(27, 840), pm(3, t)),
                ]

        self.points, self.weights = untangle(data)
        self.phi_theta = cartesian_to_spherical_sympy(self.points)
        return
