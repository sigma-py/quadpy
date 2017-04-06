# -*- coding: utf-8 -*-
#
from .helpers import _s21
import numpy


class Vertex(object):
    def __init__(self):
        self.weights = 1.0/3.0 * numpy.ones(3)
        bary = _s21(0.0)
        self.points = bary[:, [1, 2]]
        self.degree = 1
        self.name = 'vertex'
        return
