import ndim
import numpy

from .._exception import QuadpyError
from ..helpers import QuadratureScheme


class SnScheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=1.0e-14):
        self.domain = f"Sn (n={dim})"
        self.dim = dim
        assert points.shape[0] == dim
        super().__init__(name, weights, points, degree, source, tol)

    def points_inside(self):
        return numpy.einsum("ij,ij->i", self.points, self.points) < 1

    def points_inside_or_boundary(self):
        return numpy.einsum("ij,ij->i", self.points, self.points) <= 1

    def integrate(self, f, center, radius, dot=numpy.dot):
        center = numpy.array(center)
        rr = numpy.multiply.outer(radius, self.points.T)
        rr = numpy.swapaxes(rr, 0, -2)
        x = (rr + center).T
        fx = numpy.array(f(x))
        if fx.shape[-len(x.shape[1:]) :] != x.shape[1:]:
            string = ", ".join(str(val) for val in x.shape[1:])
            raise QuadpyError(
                f"Wrong return value shape {fx.shape}. " f"Expected (..., {string})."
            )
        ref_vol = ndim.nball.volume(self.dim, r=numpy.asarray(radius), symbolic=False)
        return ref_vol * dot(fx, self.weights)
