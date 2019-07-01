# -*- coding: utf-8 -*-
#
import numpy

import orthopy

from ..tools import scheme_from_rc
from ._helpers import LineSegmentScheme


def gauss_radau(n, a=0.0, b=0.0):
    assert n >= 2
    degree = 2 * n - 1
    _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.jacobi(
        n, a, b, "monic"
    )
    flt = numpy.vectorize(float)
    alpha = flt(alpha)
    beta = flt(beta)
    points, weights = _radau(alpha, beta, -1.0)
    return LineSegmentScheme("Gauss-Radau", degree, weights, points)


def _radau(alpha, beta, xr):
    """From <http://www.scientificpython.net/pyblog/radau-quadrature>:
    Compute the Radau nodes and weights with the preassigned node xr.

    Based on the section 7 of the paper

        Some modified matrix eigenvalue problems,
        Gene Golub,
        SIAM Review Vol 15, No. 2, April 1973, pp.318--334.
    """
    from scipy.linalg import solve_banded

    n = len(alpha) - 1
    f = numpy.zeros(n)
    f[-1] = beta[-1]
    A = numpy.vstack((numpy.sqrt(beta), alpha - xr))
    J = numpy.vstack((A[:, 0:-1], A[0, 1:]))
    delta = solve_banded((1, 1), J, f)
    alphar = alpha.copy()
    alphar[-1] = xr + delta[-1]
    x, w = scheme_from_rc(alphar, beta, mode="numpy")
    return x, w
