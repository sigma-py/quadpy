# -*- coding: utf-8 -*-
#
import numpy
from ..helpers import backend_to_function


class E3rScheme(object):
    def __init__(self, name, weights, points, degree, citation):
        self.name = name
        self.citation = citation
        self.degree = degree
        self.weights = weights
        self.points = points
        return

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        return dot(f(flt(self.points).T), flt(self.weights))

    def show(scheme, backend="vtk"):
        """Displays scheme for E_3^r quadrature.
        """
        backend_to_function[backend](
            scheme.points, scheme.weights, volume=8 * numpy.pi, edges=[]
        )
        return
