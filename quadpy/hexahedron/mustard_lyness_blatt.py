# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import fs_rr0, fs_r00, pm_rrr, z
from ..helpers import untangle


class MustardLynessBlatt(object):
    '''
    D. Mustard, J.N. Lyness, J.M. Blatt,
    Numerical quadrature in n dimensions,
    Comput J (1963) 6 (1): 75-87,
    <https://doi.org/10.1093/comjnl/6.1.75>.
    '''
    def __init__(self, index):
        if index == 1:
            self.degree = 3
            data = [
                (fr(1, 2), z()),
                (fr(1, 24), fs_rr0(1)),
                ]
        elif index == 2:
            self.degree = 3
            data = [
                (fr(2, 9), z()),
                (fr(1, 9), fs_r00(1)),
                (fr(1, 72), pm_rrr(1))
                ]
        elif index == 3:
            self.degree = 3
            data = [
                (+fr(1, 6), fs_rr0(1)),
                (-fr(1, 8), pm_rrr(1))
                ]
        elif index == 4:
            self.degree = 5
            data = [
                (-fr(2, 45), z()),
                (+fr(2, 45), fs_r00(1)),
                (+fr(4, 45), pm_rrr(fr(1, 2))),
                (fr(1, 120), pm_rrr(1)),
                ]
        elif index == 5:
            self.degree = 5
            data = [
                (-fr(19, 15), z()),
                (+fr(16, 45), fs_r00(fr(1, 2))),
                (-fr(1, 30), fs_r00(1)),
                (+fr(1, 36), fs_rr0(1)),
                ]
        elif index == 6:
            self.degree = 5
            data = [
                (-fr(4, 3), z()),
                (+fr(16, 45), fs_r00(fr(1, 2))),
                (fr(1, 90), fs_rr0(1)),
                (fr(1, 120), pm_rrr(1)),
                ]
        else:
            assert index == 7
            self.degree = 5
            data = [
                (fr(2, 45), z()),
                (fr(1, 45), fs_rr0(1)),
                (fr(4, 45), pm_rrr(fr(1, 2))),
                (fr(-1, 360), pm_rrr(1)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 8
        return
