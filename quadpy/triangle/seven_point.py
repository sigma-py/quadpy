# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import untangle2


class SevenPoint(object):
    def __init__(self):
        self.name = 'seven-point'
        self.degree = 3
        data = {
            's3': [fr(9, 20)],
            's2': [
                [fr(1, 20), 0],
                [fr(2, 15), fr(1, 2)],
                ],
            }
        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
