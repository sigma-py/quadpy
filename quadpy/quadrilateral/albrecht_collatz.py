# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from .helpers import _symm_r_0, _symm_s, _z, _pm, _pm2
from ..helpers import untangle


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://dx.doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self, index):
        self.name = 'AlbrechtCollatz({})'.format(index)
        if index == 1:
            self.degree = 3
            data = [
                (fr(5, 12), _z()),
                (fr(1, 8), _symm_r_0(1)),
                (fr(1, 48), _symm_s(1))
                ]
        elif index == 2:
            self.degree = 5
            r = sqrt(fr(3, 5))
            s = sqrt(fr(1, 3))
            t = sqrt(fr(14, 15))
            data = [
                (fr(5, 36), _pm2(r, s)),
                (fr(5, 63), _pm(0, t)),
                (fr(2, 7), _z())
                ]
        elif index == 3:
            self.degree = 5
            r = sqrt(fr(7, 15))
            s, t = [sqrt((7 + i*sqrt(24)) / 15) for i in [+1, -1]]
            data = [
                (fr(2, 7), _z()),
                (fr(25, 168), _pm(r, r)),
                (fr(5, 48), _pm(+s, -t)),
                (fr(5, 48), _pm(+t, -s)),
                ]
        else:
            assert index == 4
            self.degree = 5
            data = [
                (fr(2, 45), _z()),
                (fr(2, 45), _symm_r_0(1)),
                (fr(1, 60), _symm_s(1)),
                (fr(8, 45), _symm_s(fr(1, 2))),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
