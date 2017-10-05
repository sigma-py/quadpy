# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import fs_r00, fs_rr0, z
from ..helpers import untangle


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://dx.doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self):
        self.degree = 3
        data = [
            (fr(1, 4), z()),
            (fr(1, 12), fs_r00(1)),
            (fr(1, 48), fs_rr0(1)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8
        return
