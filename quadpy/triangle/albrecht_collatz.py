# -*- coding: utf-8 -*-
#
from __future__ import division

from sympy import Rational

from .helpers import _s21

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
            (Rational(2, 30), _s21(Rational(1, 2))),
            (Rational(9, 15), _s21(Rational(1, 6))),
            ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        self.weights /= 2
        return
