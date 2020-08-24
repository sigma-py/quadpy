from math import pi

import numpy

from ..helpers import QuadratureScheme, backend_to_function, expand_symmetries


class S3Scheme(QuadratureScheme):
    def __init__(self, name, source, degree, symmetry_data, tol=1.0e-14):
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data)
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = "S3"

    def show(self, backend="vtk"):
        """Displays scheme for 3D ball quadrature.
        """
        backend_to_function[backend](
            self.points,
            self.weights,
            volume=4.0 / 3.0 * pi,
            edges=[],
            balls=[((0.0, 0.0, 0.0), 1.0)],
        )
        return

    def integrate(self, f, center, radius, dot=numpy.dot):
        center = numpy.asarray(center)
        rr = numpy.multiply.outer(radius, self.points.T)
        rr = numpy.swapaxes(rr, 0, -2)
        ff = numpy.asarray(f((rr + center).T))
        return 4 / 3 * pi * numpy.asarray(radius) ** 3 * dot(ff, self.weights)
