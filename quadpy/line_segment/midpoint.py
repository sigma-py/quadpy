# -*- coding: utf-8 -*-
#
import numpy

from .helpers import LineSegmentScheme


def Midpoint():
    weights = numpy.array([2.0])
    points = numpy.array([0.0])
    return LineSegmentScheme("Midpoint rule", 1, weights, points)
