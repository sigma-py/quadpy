import numpy

from ..helpers import backend_to_function


class E3r2Scheme:
    def __init__(self, name, weights, points, degree, citation):
        self.name = name
        self.citation = citation
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
        return

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        return dot(f(flt(self.points).T), flt(self.weights))

    def show(scheme, backend="vtk"):
        """Displays scheme for E_3^{r^2} quadrature.
        """
        backend_to_function[backend](
            scheme.points, scheme.weights, volume=8 * numpy.pi, edges=[]
        )
        return
