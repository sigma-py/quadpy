# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _s21


class Vertex(object):
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y

        self.weights = numpy.full(3, frac(1, 3))
        self.bary = _s21(0)
        self.points = self.bary[:, 1:]
        self.degree = 1
        self.name = 'vertex'
        return
