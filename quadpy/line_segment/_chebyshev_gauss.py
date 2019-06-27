# -*- coding: utf-8 -*-
#
import numpy

import orthopy

from ..tools import scheme_from_rc
from ._helpers import LineSegmentScheme


def chebyshev_gauss_1(n, mode="numpy"):
    """Chebyshev-Gauss quadrature for \\int_{-1}^1 f(x) / sqrt(1+x^2) dx.
    """
    degree = n if n % 2 == 1 else n + 1

    # TODO make explicit for all modes
    if mode == "numpy":
        points = numpy.cos((2 * numpy.arange(1, n + 1) - 1) / (2 * n) * numpy.pi)
        weights = numpy.full(n, numpy.pi / n)
    else:
        _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.chebyshev1(
            n, "monic", symbolic=True
        )
        points, weights = scheme_from_rc(alpha, beta, mode)
    return LineSegmentScheme("Chebyshev-Gauss 1", degree, weights, points)


def chebyshev_gauss_2(n, mode="numpy", decimal_places=None):
    """Chebyshev-Gauss quadrature for \\int_{-1}^1 f(x) * sqrt(1+x^2) dx.
    """
    degree = n if n % 2 == 1 else n + 1

    # TODO make explicit for all modes
    if mode == "numpy":
        points = numpy.cos(numpy.pi * numpy.arange(1, n + 1) / (n + 1))
        weights = (
            numpy.pi
            / (n + 1)
            * (numpy.sin(numpy.pi * numpy.arange(1, n + 1) / (n + 1))) ** 2
        )
    else:
        _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.chebyshev2(
            n, "monic", symbolic=True
        )
        points, weights = scheme_from_rc(alpha, beta, mode)
    return LineSegmentScheme("Chebyshev-Gauss 2", degree, weights, points)
