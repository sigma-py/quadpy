import ndim
import numpy

from ..helpers import QuadratureScheme


class EnrScheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = f"Enr (n={dim})"
        self.dim = dim

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = ndim.enr.volume(self.dim)
        return ref_vol * dot(f(flt(self.points).T), flt(self.weights))
