from math import pi, sqrt

import numpy

from ..helpers import QuadratureScheme, backend_to_function


class E3r2Scheme(QuadratureScheme):
    def __init__(self, name, weights, points, degree, source, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = "E3r2"

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = sqrt(pi) ** 3
        return ref_vol * dot(f(flt(self.points)), flt(self.weights))

    def show(scheme, backend="vtk"):
        """Displays scheme for E_3^{r^2} quadrature.
        """
        backend_to_function[backend](
            scheme.points, scheme.weights, volume=8 * numpy.pi, edges=[]
        )
