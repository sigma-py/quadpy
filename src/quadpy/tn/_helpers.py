from __future__ import annotations

import math
from typing import Callable

import numpy as np

from .._exception import QuadpyError
from ..helpers import QuadratureScheme


class TnScheme(QuadratureScheme):
    def __init__(
        self,
        name: str,
        dim: int,
        weights,
        points,
        degree: int,
        source,
        tol: float = 1.0e-14,
        comments: list[str] | None = None,
    ):
        self.domain = f"Tn (n={dim})"
        self.dim = dim
        super().__init__(name, weights, points, degree, source, tol, comments)

    def points_inside(self):
        return np.all((0 < self.points) & (self.points < 1))

    def points_inside_or_boundary(self):
        return np.all((0 <= self.points) & (self.points <= 1))

    def integrate(self, f: Callable, simplex, dot=np.dot):
        flt = np.vectorize(float)
        simplex = np.asarray(simplex)
        if simplex.shape[0] != self.dim + 1 or simplex.shape[-1] < self.dim:
            string = ", ".join(str(val) for val in simplex.shape)
            raise QuadpyError(
                f"Wrong domain shape. Expected ({self.dim + 1}, ..., n >= {self.dim}), "
                f"got ({string})."
            )

        x = transform(flt(self.points), simplex.T)
        vol = get_vol(simplex)

        fx = np.asarray(f(x))
        if fx.shape[-len(x.shape[1:]) :] != x.shape[1:]:
            string = ", ".join(str(val) for val in x.shape[1:])
            raise QuadpyError(
                f"Wrong return value shape {fx.shape}. Expected (..., {string})."
            )

        return vol * dot(fx, flt(self.weights))


def transform(points, simplex):
    """Transform the points `xi` from the reference simplex onto `simplex`."""
    # For n == 2:
    # x = (
    #     + outer(triangle[0].T, 1.0 - xi[0] - xi[1])
    #     + outer(triangle[1].T, xi[0])
    #     + outer(triangle[2].T, xi[1])
    #     )
    return np.dot(simplex, points)


def get_vol(simplex):
    # Compute the volume via the Cayley-Menger determinant
    # <http://mathworld.wolfram.com/Cayley-MengerDeterminant.html>. One advantage is
    # that it can compute the volume of the simplex indenpendent of the dimension of the
    # space in which it is embedded.

    # compute all edge lengths
    edges = np.subtract(simplex[:, None], simplex[None, :])
    ei_dot_ej = np.einsum("...k,...k->...", edges, edges)

    j = simplex.shape[0] - 1
    a = np.empty((j + 2, j + 2) + ei_dot_ej.shape[2:])
    a[1:, 1:] = ei_dot_ej
    a[0, 1:] = 1.0
    a[1:, 0] = 1.0
    a[0, 0] = 0.0

    a = np.moveaxis(a, (0, 1), (-2, -1))
    det = np.linalg.det(a)

    vol = np.sqrt((-1.0) ** (j + 1) / 2 ** j / math.factorial(j) ** 2 * det)
    return vol
