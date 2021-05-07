import ndim
import numpy as np

from ..helpers import QuadratureScheme


class UnScheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = f"Un (n={dim})"
        self.dim = dim

    def integrate(self, f, center, radius, dot=np.dot):
        center = np.array(center)
        rr = np.multiply.outer(radius, self.points)
        rr = np.swapaxes(rr, 0, -2)
        ff = np.array(f((rr + center).T))
        ref_vol = ndim.nsphere.volume(self.dim, r=radius)
        return ref_vol * dot(ff, self.weights)
