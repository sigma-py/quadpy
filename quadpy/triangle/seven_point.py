# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import untangle2


class SevenPoint(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.name = "seven-point"
        self.degree = 3
        data = {
            "s3": [frac(9, 20)],
            "s2": [[frac(1, 20), 0], [frac(2, 15), frac(1, 2)]],
        }
        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return
