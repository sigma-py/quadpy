import itertools
from typing import Callable

import numpy as np

from .._exception import QuadpyError
from ..helpers import QuadratureScheme, n_outer


class CnScheme(QuadratureScheme):
    def __init__(
        self,
        name: str,
        dim: int,
        weights,
        points,
        degree: int,
        source,
        tol=1.0e-14,
        comments=None,
    ):
        self.domain = f"Cn (n={dim})"
        self.dim = dim
        assert points.shape[0] == dim, f"points.shape == {points.shape}, dim = {dim}"
        super().__init__(name, weights, points, degree, source, tol, comments)

    def integrate(self, f: Callable, ncube, dot=np.dot):
        ncube = np.asarray(ncube)
        if (
            ncube.shape[: self.dim] != tuple(self.dim * [2])
            or ncube.shape[-1] < self.dim
        ):
            expected = ", ".join(str(i) for i in self.dim * [2])
            string = ", ".join(str(val) for val in ncube.shape)
            raise QuadpyError(
                f"Wrong domain shape. Expected ({expected}, ..., n >= {self.dim}), "
                f"got ({string})."
            )
        x = transform(self.points, ncube).T
        detJ = get_detJ(self.points, ncube)

        fx = np.asarray(f(x))
        if fx.shape[-len(x.shape[1:]) :] != x.shape[1:]:
            string = ", ".join(str(val) for val in x.shape[1:])
            raise QuadpyError(
                f"Wrong return value shape {fx.shape}. Expected (..., {string})."
            )

        ref_vol = 2 ** np.prod(self.dim)
        return ref_vol * dot(fx * abs(detJ), self.weights)

    def points_inside(self):
        """Are all points strictly inside the domain?"""
        return np.all((-1 < self.points) & (self.points < 1))

    def points_inside_or_boundary(self):
        """Are all points inside the domain or on its boundary?"""
        return np.all((-1 <= self.points) & (self.points <= 1))


def transform(xi, cube):
    """Transform the points `xi` from the reference cube to `cube`."""
    # For d==2, the result used to be computed with
    #
    # out = (
    #     + outer(0.25*(1.0-xi[0])*(1.0-xi[1]), cube[0, 0])
    #     + outer(0.25*(1.0+xi[0])*(1.0-xi[1]), cube[1, 0])
    #     + outer(0.25*(1.0-xi[0])*(1.0+xi[1]), cube[0, 1])
    #     + outer(0.25*(1.0+xi[0])*(1.0+xi[1]), cube[1, 1])
    #     )
    #
    # This array of multiplications and additions is reminiscent of dot(), and
    # indeed tensordot() can handle the situation. We just need to compute the
    # `1+-xi` products and align them with `cube`.
    one_mp_xi = np.stack([0.5 * (1.0 - xi), 0.5 * (1.0 + xi)], axis=1)
    a = n_outer(one_mp_xi)

    # <https://stackoverflow.com/q/45372098/353337>
    d = xi.shape[0]
    return np.tensordot(a, cube, axes=(range(d), range(d)))


def get_detJ(xi, cube):
    """Get the determinant of the transformation matrix."""
    # For d==2, the result can be computed with
    # ```
    # J0 = (
    #     - np.multiply.outer(0.25*(1-xi[1]), quad[0, 0])
    #     + np.multiply.outer(0.25*(1-xi[1]), quad[1, 0])
    #     - np.multiply.outer(0.25*(1+xi[1]), quad[0, 1])
    #     + np.multiply.outer(0.25*(1+xi[1]), quad[1, 1])
    #     ).T
    # J1 = (
    #     - np.multiply.outer(0.25*(1-xi[0]), quad[0, 0])
    #     - np.multiply.outer(0.25*(1+xi[0]), quad[1, 0])
    #     + np.multiply.outer(0.25*(1-xi[0]), quad[0, 1])
    #     + np.multiply.outer(0.25*(1+xi[0]), quad[1, 1])
    #     ).T
    # out = J0[0]*J1[1] - J1[0]*J0[1]
    # ```
    # Like transform(), simplify here and form the determinant explicitly.
    d = xi.shape[0]

    one_mp_xi = np.stack([0.5 * (1.0 - xi), 0.5 * (1.0 + xi)], axis=1)

    # Build the Jacobi matrix row by row.
    J = []
    for k in range(d):
        a = one_mp_xi.copy()
        a[k, 0, :] = -0.5
        a[k, 1, :] = +0.5
        a0 = n_outer(a)
        J.append(np.tensordot(a0, cube, axes=(range(d), range(d))).T)

    # `det` needs the square at the end. Fortran...
    # For d==2 or d==3, we could avoid this copy and compute the determinant
    # with their elementary formulas, i.e.,
    #
    #     + J[0][0]*J[1][1] - J[1][0]*J[0][1];
    #
    #     + J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2]
    #     - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0].
    #
    J = np.array(J)
    J = np.moveaxis(J, (0, 1), (-2, -1))
    out = np.linalg.det(J)
    return out


def integrate_monomial_over_ncube(ncube_limits, exp):
    return np.prod(
        [
            (a[1] ** (k + 1) - a[0] ** (k + 1)) / (k + 1)
            for a, k in zip(ncube_limits, exp)
        ]
    )


def ncube_points(*xyz):
    """Given the end points of an n-cube aligned with the coordinate axes, this returns
    the corner points of the cube in the correct data structure.
    """
    return np.moveaxis(np.array(np.meshgrid(*xyz, indexing="ij")), 0, -1)


def _fs11(n, r, s):
    """Get all permutations of [+-r, +-s, ..., +-s] of length n.
    len(out) == n * 2**n.
    """
    # <https://stackoverflow.com/a/45321972/353337>
    k1 = 1
    k2 = n - 1
    idx = itertools.combinations(range(k1 + k2), k1)
    vs = ((s if j not in i else r for j in range(k1 + k2)) for i in idx)
    return np.array(
        list(
            itertools.chain.from_iterable(
                itertools.product(*((+vij, -vij) for vij in vi)) for vi in vs
            )
        )
    )


def _s(n, a, b):
    """Get all permutations of [a, b, ..., b] of length n. len(out) == n."""
    out = np.full((n, n), b)
    np.fill_diagonal(out, a)
    return out


def _s2(n, a):
    """Get all permutations of [a, a, 0, 0, ..., 0] of length n.
    len(out) == (n-1)*n / 2.
    """
    return _s11(n, a, a)


def _s11(n, a, b):
    """Get all permutations of [a, b, 0, 0, ..., 0] of length n.
    len(out) == (n-1)*n.
    """
    s = [a, b] + (n - 2) * [0]
    # First permutations, then set can be really inefficient if items are
    # repeated. Check out <https://stackoverflow.com/q/6284396/353337> for
    # improvements.
    return np.array(list(set(itertools.permutations(s, n))))
