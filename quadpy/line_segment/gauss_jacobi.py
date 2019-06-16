# -*- coding: utf-8 -*-
#
import orthopy

from ..tools import scheme_from_rc
from .helpers import LineSegmentScheme


def GaussJacobi(n, alpha, beta, mode="numpy", decimal_places=None):
    """
    Gauss-Jacobi quadrature.
    """
    degree = 2 * n - 1

    _, _, a, b = orthopy.line_segment.recurrence_coefficients.jacobi(
        n, alpha, beta, "monic", symbolic=True
    )
    points, weights = scheme_from_rc(
        a, b, mode=mode, decimal_places=decimal_places
    )
    return LineSegmentScheme("Gauss-Jacobi", degree, weights, points)
