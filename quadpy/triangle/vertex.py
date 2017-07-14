# -*- coding: utf-8 -*-
#
import numpy
from .helpers import _s21


class Vertex(object):
    def __init__(self):
        self.weights = numpy.full(3, 1.0/3.0)
        bary = _s21(0.0)
        self.points = bary[:, [1, 2]]
        self.degree = 1
        self.name = 'vertex'
        return
