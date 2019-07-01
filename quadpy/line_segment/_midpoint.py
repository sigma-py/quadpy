import numpy

from ._helpers import LineSegmentScheme


def midpoint():
    weights = numpy.array([2.0])
    points = numpy.array([0.0])
    return LineSegmentScheme("Midpoint rule", 1, weights, points)
