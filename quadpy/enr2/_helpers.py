import ndim
import numpy

from ..helpers import QuadratureScheme


class Enr2Scheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=2.7e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = f"Enr2 (n={dim}, {tol})"
        self.dim = dim

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = ndim.enr2.volume(self.dim, "physicists")
        return ref_vol * dot(f(flt(self.points).T), flt(self.weights))
