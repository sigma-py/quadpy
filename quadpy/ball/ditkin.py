# -*- coding: utf-8 -*-
#
from sympy import Rational as fr, sqrt, pi

from ..helpers import untangle, z, pm_array0, pm


class Ditkin(object):
    '''
    V. A. Ditkin,
    On certain approximate formulas for the calculation of triple integrals,
    Doklady Akad. Nauk SSSR (N.S.) 62 (1948), 445–447 (Russian).
    '''
    def __init__(self, index, alpha=0):
        if index == 1:
            self.degree = 5

            B0 = fr(4, (alpha + 5)**2)
            B1 = fr((alpha + 3) * (alpha + 7), 12 * (alpha + 5)**2)

            r, s = [
                sqrt((alpha+5) * (5+t*sqrt(5)) / 10 / (alpha+7))
                for t in [+1, -1]
                ]

            data = [
                (B0, z(3)),
                (B1, pm_array0(3, [r, s], [0, 1])),
                (B1, pm_array0(3, [r, s], [1, 2])),
                (B1, pm_array0(3, [r, s], [2, 0])),
                ]
        elif index == 2:
            self.degree = 5

            B0 = fr(4, 25)
            B1 = fr(21, 500)

            r, s = [sqrt((15 + 5*t*sqrt(5)) / 42) for t in [+1, -1]]
            t = sqrt(fr(5, 21))

            data = [
                (B0, z(3)),
                (B1, pm_array0(3, [r, s], [0, 1])),
                (B1, pm_array0(3, [r, s], [1, 2])),
                (B1, pm_array0(3, [r, s], [2, 0])),
                (B1, pm(3, t)),
                ]
        else:
            assert index == 3, 'Illegal Ditkin index {}.'.format(index)
            self.degree = 7

            B0 = fr(16, 175)
            B1 = fr(81, 1400)
            B2 = fr(3, 280)

            sqrt5 = sqrt(5)
            r, s = [sqrt((5 + pm_*sqrt5) / 18) for pm_ in [+1, -1]]
            t = sqrt(fr(1, 3))
            u, v = [sqrt((3 - pm_*sqrt5) / 6) for pm_ in [+1, -1]]

            data = [
                (B0, z(3)),
                (B1, pm_array0(3, [r, s], [0, 1])),
                (B1, pm_array0(3, [r, s], [1, 2])),
                (B1, pm_array0(3, [r, s], [2, 0])),
                (B2, pm_array0(3, [u, v], [0, 1])),
                (B2, pm_array0(3, [u, v], [1, 2])),
                (B2, pm_array0(3, [u, v], [2, 0])),
                (B2, pm(3, t)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= fr(4, 3) * pi
        return
