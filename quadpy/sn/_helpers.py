import ndim
import numpy

from ..helpers import QuadratureScheme


class SnScheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = f"Sn (n={dim})"
        self.dim = dim

    def points_inside(self):
        return numpy.einsum("ij,ij->i", self.points, self.points) < 1

    def points_inside_or_boundary(self):
        return numpy.einsum("ij,ij->i", self.points, self.points) <= 1

    def integrate(self, f, center, radius, dot=numpy.dot):
        center = numpy.array(center)
        rr = numpy.multiply.outer(radius, self.points)
        rr = numpy.swapaxes(rr, 0, -2)
        ff = numpy.array(f((rr + center).T))
        ref_vol = ndim.nball.volume(self.dim, r=numpy.asarray(radius), symbolic=False)
        return ref_vol * dot(ff, self.weights)
