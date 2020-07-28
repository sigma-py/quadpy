import ndim
import numpy

from ..helpers import QuadratureScheme


class EnrScheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=1.0e-14):
        self.domain = f"Enr (n={dim})"
        self.dim = dim
        super().__init__(name, weights, points, degree, source, tol)

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = ndim.enr.volume(self.dim)
        return ref_vol * dot(f(flt(self.points).T), flt(self.weights))
