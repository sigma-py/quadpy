# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from .helpers import fs_r00, fs_rr0, pm_rrr
from ..helpers import untangle


class HammerWymore(object):
    '''
    Preston C. Hammer and A. Wayne Wymore,
    Numerical evaluation of multiple integrals. I,
    Math. Comp. 11 (1957), 59-67,
    <https://doi.org/10.1090/S0025-5718-1957-0087220-6>.
    '''
    def __init__(self):
        self.degree = 7

        r2 = fr(6, 7)
        s2 = (960 - 3*sqrt(28798)) / 2726
        t2 = (960 + 3*sqrt(28798)) / 2726

        r = sqrt(r2)
        s = sqrt(s2)
        t = sqrt(t2)

        B1 = fr(1078, 29160)
        B2 = fr(343, 29160)
        B3 = (774*t2 - 230) / (9720 * (t2-s2))
        B4 = (230 - 774*s2) / (9720 * (t2-s2))

        data = [
            (B1, fs_r00(r)),
            (B2, fs_rr0(r)),
            (B3, pm_rrr(s)),
            (B4, pm_rrr(t)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8
        return
