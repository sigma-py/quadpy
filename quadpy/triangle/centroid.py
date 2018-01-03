# -*- coding: utf-8 -*-
#
import numpy
from .helpers import _s3


class Centroid(object):
    def __init__(self, symbolic=False):
        self.weights = numpy.array([1])
        self.bary = _s3(symbolic)
        self.points = self.bary[:, 1:]
        self.degree = 1
        self.name = 'centroid'
        return
