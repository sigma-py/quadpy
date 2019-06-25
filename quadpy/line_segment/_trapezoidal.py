# -*- coding: utf-8 -*-
#
import numpy

from ._helpers import LineSegmentScheme


def trapezoidal():
    weights = numpy.array([1.0, 1.0])
    points = numpy.array([-1.0, 1.0])
    return LineSegmentScheme("Trapezoidal rule", 1, weights, points)
