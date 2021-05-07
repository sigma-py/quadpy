from typing import Callable

import ndim
import numpy as np

from .._exception import QuadpyError
from ..helpers import QuadratureScheme


class EnrScheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=1.0e-14):
        self.domain = f"Enr (n={dim})"
        self.dim = dim
        assert points.shape[0] == dim
        super().__init__(name, weights, points, degree, source, tol)

    def integrate(self, f: Callable, dot=np.dot):
        flt = np.vectorize(float)
        ref_vol = ndim.enr.volume(self.dim)

        x = flt(self.points)
        fx = np.asarray(f(x))

        assert len(x.shape) == 2
        if len(x.shape) == 2 and fx.shape[-1] != x.shape[1]:
            raise QuadpyError(
                f"Wrong return value shape {fx.shape}. Expected (:, {x.shape[1]})."
            )

        return ref_vol * dot(fx, flt(self.weights))
