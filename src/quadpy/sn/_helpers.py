import ndim
import numpy as np

from .._exception import QuadpyError
from ..helpers import QuadratureScheme


class SnScheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=1.0e-14):
        self.domain = f"Sn (n={dim})"
        self.dim = dim
        assert points.shape[0] == dim
        super().__init__(name, weights, points, degree, source, tol)

    def points_inside(self):
        return np.einsum("ij,ij->i", self.points, self.points) < 1

    def points_inside_or_boundary(self):
        return np.einsum("ij,ij->i", self.points, self.points) <= 1

    def integrate(self, f, center, radius, dot=np.dot):
        center = np.array(center)
        rr = np.multiply.outer(radius, self.points.T)
        rr = np.swapaxes(rr, 0, -2)
        x = (rr + center).T
        fx = np.array(f(x))
        if fx.shape[-len(x.shape[1:]) :] != x.shape[1:]:
            string = ", ".join(str(val) for val in x.shape[1:])
            raise QuadpyError(
                f"Wrong return value shape {fx.shape}. Expected (..., {string})."
            )
        ref_vol = ndim.nball.volume(self.dim, r=np.asarray(radius), symbolic=False)
        return ref_vol * dot(fx, self.weights)
