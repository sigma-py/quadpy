# -*- coding: utf-8 -*-
#
import numpy
import orthopy

from ..tools import scheme_from_rc
from .helpers import LineSegmentScheme


def GaussLegendre(n, mode="numpy", decimal_places=None):
    """
    Gauss-Legendre quadrature.
    """
    degree = 2 * n - 1

    if mode == "numpy":
        points, weights = numpy.polynomial.legendre.leggauss(n)
    else:
        _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.legendre(
            n, "monic", symbolic=True
        )
        points, weights = scheme_from_rc(
            alpha, beta, mode=mode, decimal_places=decimal_places
        )
    return LineSegmentScheme("Gauss-Legendre", degree, weights, points)
