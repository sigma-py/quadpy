import numpy
import orthopy

from ..tools import scheme_from_rc
from ._helpers import C1Scheme


def gauss_legendre(n, mode="numpy"):
    degree = 2 * n - 1

    if mode == "numpy":
        points, weights = numpy.polynomial.legendre.leggauss(n)
    else:
        rc = orthopy.c1.legendre.RecurrenceCoefficients("monic", symbolic=True)
        _, alpha, beta = numpy.array([rc[k] for k in range(n)]).T
        beta[0] = rc.int_1
        points, weights = scheme_from_rc(alpha, beta, mode=mode)
    return C1Scheme("Gauss-Legendre", degree, weights, points)
