# -*- coding: utf-8 -*-
#
from math import pi
import numpy

from .. import helpers


# @dataclass
class BallScheme(object):
    # name: str
    # degree: int
    # weights: list
    # points: list

    def __init__(self, name, citation, degree, weights, points):
        self.name = name
        self.citation = citation
        self.degree = degree
        self.weights = weights
        self.points = points
        return

    def show(self, backend="vtk"):
        """Displays scheme for 3D ball quadrature.
        """
        helpers.backend_to_function[backend](
            self.points,
            self.weights,
            volume=4.0 / 3.0 * pi,
            edges=[],
            balls=[((0.0, 0.0, 0.0), 1.0)],
        )
        return

    def integrate(self, f, center, radius, dot=numpy.dot):
        center = numpy.asarray(center)
        rr = numpy.multiply.outer(radius, self.points)
        rr = numpy.swapaxes(rr, 0, -2)
        ff = numpy.asarray(f((rr + center).T))
        return numpy.asarray(radius) ** 3 * dot(ff, self.weights)
