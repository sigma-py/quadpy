import numpy as np
import orthopy

from ..tools import scheme_from_rc
from ._helpers import C1Scheme


def gauss_legendre(n, mode="numpy"):
    degree = 2 * n - 1

    if mode == "numpy":
        points, weights = np.polynomial.legendre.leggauss(n)
    else:
        rc = orthopy.c1.legendre.RecurrenceCoefficients("monic", symbolic=True)
        _, alpha, beta = np.array([rc[k] for k in range(n)]).T
        points, weights = scheme_from_rc(alpha, beta, rc.int_1, mode=mode)
    return C1Scheme("Gauss-Legendre", degree, weights, points)
