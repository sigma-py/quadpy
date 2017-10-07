# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from .helpers import _symm_r_0, _symm_s, _z
from ..helpers import untangle


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://dx.doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self, index):
        self.name = 'Tyler({})'.format(index)
        if index == 1:
            self.degree = 5
            data = [
                (-fr(28, 45), _z()),
                (fr(1, 36), _symm_s(1)),
                (fr(1, 45), _symm_r_0(1)),
                (fr(16, 45), _symm_r_0(fr(1, 2))),
                ]
        elif index == 2:
            self.degree = 7
            r = sqrt(fr(6, 7))
            s, t = [sqrt((114 - i*3*sqrt(583)) / 287) for i in [+1, -1]]
            B1 = fr(49, 810)
            B2, B3 = [
                (178981 + i * 2769 * sqrt(583)) / 1888920
                for i in [+1, -1]
                ]
            data = [
                (B1, _symm_r_0(r)),
                (B2, _symm_s(s)),
                (B3, _symm_s(t)),
                ]
        else:
            assert index == 3
            self.degree = 7
            r = fr(2, 3)
            s = fr(1, 3)
            t = fr(1, 2)
            data = [
                (fr(449, 315), _z()),
                (fr(37, 1260), _symm_r_0(1)),
                (fr(3, 28), _symm_r_0(r)),
                (-fr(69, 140), _symm_r_0(s)),
                (fr(7, 540), _symm_s(1)),
                (fr(32, 135), _symm_s(t)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
