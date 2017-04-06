# -*- coding: utf-8 -*-
#
from .helpers import _s3
import numpy


class Centroid(object):
    def __init__(self):
        self.weights = numpy.array([1.0])
        bary = _s3()
        self.points = bary[:, [1, 2]]
        self.degree = 1
        self.name = 'centroid'
        return
