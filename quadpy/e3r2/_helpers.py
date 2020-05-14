from math import pi, sqrt

import numpy

from ..helpers import QuadratureScheme, backend_to_function


class E3r2Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree, source):
        self.domain = "E3r2"
        self.name = name
        self.source = source
        self.degree = degree

        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int_]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int_]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = sqrt(pi) ** 3
        return ref_vol * dot(f(flt(self.points).T), flt(self.weights))

    def show(scheme, backend="vtk"):
        """Displays scheme for E_3^{r^2} quadrature.
        """
        backend_to_function[backend](
            scheme.points, scheme.weights, volume=8 * numpy.pi, edges=[]
        )
