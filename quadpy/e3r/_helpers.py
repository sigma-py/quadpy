from math import pi

import numpy

from ..helpers import QuadratureScheme, backend_to_function


class E3rScheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree, source, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = "E3r"

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = 8 * pi
        return ref_vol * dot(f(flt(self.points).T), flt(self.weights))

    def show(scheme, backend="vtk"):
        """Displays scheme for E_3^r quadrature.
        """
        backend_to_function[backend](
            scheme.points, scheme.weights, volume=8 * numpy.pi, edges=[]
        )
