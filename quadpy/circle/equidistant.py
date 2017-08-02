# -*- coding: utf-8 -*-
#
import numpy


class Equidistant(object):
    def __init__(self, n):
        self.weights = numpy.full(n, 2 * numpy.pi / n)
        self.points = numpy.column_stack([
            numpy.cos(2*numpy.pi * numpy.arange(n) / n),
            numpy.sin(2*numpy.pi * numpy.arange(n) / n),
            ])
        self.degree = n - 1
        return
