# -*- coding: utf-8 -*-
#
import numpy
from sympy import Rational as fr

from .helpers import _s21


class Vertex(object):
    def __init__(self):
        self.weights = numpy.full(3, fr(1, 3))
        self.bary = _s21(0)
        self.points = self.bary[:, 1:]
        self.degree = 1
        self.name = 'vertex'
        return
