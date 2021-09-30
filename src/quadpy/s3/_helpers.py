from __future__ import annotations

from math import pi
from typing import Callable

import numpy as np

from .._exception import QuadpyError
from ..helpers import QuadratureScheme, backend_to_function, expand_symmetries

schemes = {}


def register(in_schemes):
    for scheme in in_schemes:
        schemes[scheme.__name__] = scheme


class S3Scheme(QuadratureScheme):
    def __init__(self, name, symmetry_data, degree, source=None, tol=1.0e-14):
        self.symmetry_data = symmetry_data
        points, weights = expand_symmetries(symmetry_data, dim=3)
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = "S3"

    def show(self, backend="vtk", **kwargs):
        """Displays scheme for 3D ball quadrature."""
        return backend_to_function[backend](
            self.points,
            self.weights,
            volume=4.0 / 3.0 * pi,
            edges=[],
            balls=[((0.0, 0.0, 0.0), 1.0)],
            **kwargs,
        )

    def integrate(self, f: Callable, center, radius, dot=np.dot):
        center = np.asarray(center)
        rr = np.multiply.outer(radius, self.points.T)
        rr = np.swapaxes(rr, 0, -2)
        x = (rr + center).T
        ff = np.asarray(f(x))
        if ff.shape[-len(x.shape[1:]) :] != x.shape[1:]:
            string = ", ".join(str(val) for val in x.shape[1:])
            raise QuadpyError(
                f"Wrong return value shape {ff.shape}. Expected (..., {string})."
            )
        return 4 / 3 * pi * np.asarray(radius) ** 3 * dot(ff, self.weights)


def get_good_scheme(degree: int) -> S3Scheme | None:
    if degree <= 7:
        return {
            0: schemes["midpoint"],
            1: schemes["midpoint"],
            2: schemes["hammer_stroud_11_3"],
            3: schemes["hammer_stroud_11_3"],
            4: schemes["hammer_stroud_12_3"],
            5: schemes["hammer_stroud_12_3"],
            6: schemes["mysovskih"],
            7: schemes["mysovskih"],
            8: schemes["stroud_14_1"],
            9: schemes["stroud_14_1"],
            10: schemes["stroud_14_1"],
            11: schemes["stroud_14_1"],
            12: schemes["stroud_14_1"],
            13: schemes["stroud_14_1"],
            14: schemes["stroud_14_1"],
        }[degree]()
    return None
