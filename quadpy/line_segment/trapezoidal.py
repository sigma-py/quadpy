# -*- coding: utf-8 -*-
#
import numpy


class Trapezoidal(object):
    def __init__(self):
        self.weights = numpy.array([1.0, 1.0])
        self.points = numpy.array([-1.0, 1.0])
        self.degree = 1
        return
