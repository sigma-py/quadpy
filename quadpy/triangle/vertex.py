# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import untangle2


class Vertex(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        data = {"s2": [[frac(1, 3), 0]]}

        self.bary, self.weights = untangle2(data)

        self.points = self.bary[:, 1:]
        self.degree = 1
        self.name = "vertex"
        return
