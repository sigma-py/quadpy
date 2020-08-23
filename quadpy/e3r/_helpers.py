from math import pi

import numpy

from ..helpers import QuadratureScheme, backend_to_function, expand_symmetries

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class E3rScheme(QuadratureScheme):
    def __init__(self, name, symmetry_data, degree, source, tol=1.0e-14):
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data)
        self.domain = "E3r"
        super().__init__(name, weights, points, degree, source, tol)

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = 8 * pi
        return ref_vol * dot(f(flt(self.points)), flt(self.weights))

    def show(scheme, backend="vtk"):
        """Displays scheme for E_3^r quadrature.
        """
        backend_to_function[backend](
            scheme.points, scheme.weights, volume=8 * numpy.pi, edges=[]
        )
