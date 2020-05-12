import numpy

from ._helpers import C1Scheme


def midpoint():
    weights = numpy.array([2.0])
    points = numpy.array([0.0])
    return C1Scheme("Midpoint rule", 1, weights, points)
