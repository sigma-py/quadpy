import math

import numpy
import scipy.special
import sympy

from ..helpers import prod


class NSimplexScheme:
    def __init__(self, name, dim, weights, points, degree, citation):
        self.name = name
        self.dim = dim
        self.degree = degree
        self.citation = citation

        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int_]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int_]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points
        return

    def integrate(self, f, simplex, dot=numpy.dot):
        flt = numpy.vectorize(float)
        simplex = numpy.asarray(simplex)
        print(self.points.shape)
        x = transform(flt(self.points).T, simplex.T)
        vol = get_vol(simplex)

        fx = numpy.asarray(f(x))

        assert (
            x.shape[1:] == fx.shape[-len(x.shape[1:]) :]
        ), "Illegal shape of f(x) (expected (..., {}), got {})".format(
            ", ".join([str(k) for k in x.shape[1:]]), fx.shape
        )
        return vol * dot(fx, flt(self.weights))


def transform(points, simplex):
    """Transform the points `xi` from the reference simplex onto `simplex`.
    """
    # For n == 2:
    # x = (
    #     + outer(triangle[0].T, 1.0 - xi[0] - xi[1])
    #     + outer(triangle[1].T, xi[0])
    #     + outer(triangle[2].T, xi[1])
    #     )
    return numpy.dot(simplex, points)


def get_vol(simplex):
    # Compute the volume via the Cayley-Menger determinant
    # <http://mathworld.wolfram.com/Cayley-MengerDeterminant.html>. One advantage is
    # that it can compute the volume of the simplex indenpendent of the dimension of the
    # space in which it is embedded.

    # compute all edge lengths
    edges = numpy.subtract(simplex[:, None], simplex[None, :])
    ei_dot_ej = numpy.einsum("...k,...k->...", edges, edges)

    j = simplex.shape[0] - 1
    a = numpy.empty((j + 2, j + 2) + ei_dot_ej.shape[2:])
    a[1:, 1:] = ei_dot_ej
    a[0, 1:] = 1.0
    a[1:, 0] = 1.0
    a[0, 0] = 0.0

    a = numpy.moveaxis(a, (0, 1), (-2, -1))
    det = numpy.linalg.det(a)

    vol = numpy.sqrt((-1.0) ** (j + 1) / 2 ** j / math.factorial(j) ** 2 * det)
    return vol


def integrate_monomial_over_unit_simplex(k, symbolic=False):
    """The integrals of monomials over the standard triangle and tetrahedron are
    given by

    \\int_T x_0^k0 * x1^k1 = (k0!*k1!) / (2+k0+k1)!,
    \\int_T x_0^k0 * x1^k1 * x2^k2 = (k0!*k1!*k2!) / (4+k0+k1+k2)!,

    see, e.g.,
    A set of symmetric quadrature rules on triangles and tetrahedra,
    Linbo Zhang, Tao Cui and Hui Liu,
    Journal of Computational Mathematics,
    Vol. 27, No. 1 (January 2009), pp. 89-96,
    <https://www.jstor.org/stable/43693493>.

    See, e.g., <https://math.stackexchange.com/q/207073/36678> for a formula in
    all dimensions.
    """
    if symbolic:
        return sympy.Rational(
            prod([math.factorial(kk) for kk in k]), math.factorial(sum(k) + len(k))
        )
    # exp-log to account for large values in numerator and denominator
    # import scipy.special
    return math.exp(
        math.fsum([scipy.special.gammaln(kk + 1) for kk in k])
        - scipy.special.gammaln(sum([kk + 1 for kk in k]) + 1)
    )
