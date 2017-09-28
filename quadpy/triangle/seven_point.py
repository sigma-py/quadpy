# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import _s3, _s21
from ..helpers import untangle


class SevenPoint(object):
    def __init__(self):
        self.name = 'seven-point'
        self.degree = 3
        data = [
            (fr(9, 20), _s3()),
            (fr(1, 20), _s21(0)),
            (fr(2, 15), _s21(fr(1, 2))),
            ]
        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
