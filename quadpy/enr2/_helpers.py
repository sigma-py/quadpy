import ndim
import numpy as np

from ..helpers import QuadratureScheme


class Enr2Scheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=2.7e-14):
        self.domain = f"Enr2 (n={dim}, {tol})"
        self.dim = dim
        assert points.shape[0] == dim
        super().__init__(name, weights, points, degree, source, tol)

    def integrate(self, f, dot=np.dot):
        flt = np.vectorize(float)
        ref_vol = ndim.enr2.volume(self.dim, "physicists")
        return ref_vol * dot(f(flt(self.points)), flt(self.weights))
