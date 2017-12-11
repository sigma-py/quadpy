# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import untangle2


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self):
        self.degree = 3

        self.data = {
            's2': [
                [fr(2, 30), fr(1, 2)],
                [fr(9, 15), fr(1, 6)],
                ],
            }
        self.bary, self.weights = untangle2(self.data)

        self.points = self.bary[:, 1:]
        self.weights /= 2
        return
