# -*- coding: utf-8 -*-
#
import numpy


class Midpoint(object):
    def __init__(self):
        self.weights = numpy.array([2.0])
        self.points = numpy.array([
            0.0
            ])
        self.degree = 1
        return
