import numpy

import orthopy

from ..tools import scheme_from_rc
from ._helpers import LineSegmentScheme


def gauss_legendre(n, mode="numpy"):
    degree = 2 * n - 1

    if mode == "numpy":
        points, weights = numpy.polynomial.legendre.leggauss(n)
    else:
        _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.legendre(
            n, "monic", symbolic=True
        )
        points, weights = scheme_from_rc(alpha, beta, mode=mode)
    return LineSegmentScheme("Gauss-Legendre", degree, weights, points)
