import numpy

from ._helpers import C1Scheme


def trapezoidal():
    weights = numpy.array([1.0, 1.0])
    points = numpy.array([-1.0, 1.0])
    return C1Scheme("Trapezoidal rule", 1, weights, points)
