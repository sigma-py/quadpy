import numpy as np
import orthopy

from ..tools import scheme_from_rc
from ._helpers import C1Scheme


def gauss_jacobi(n, alpha, beta, mode="numpy"):
    degree = 2 * n - 1

    rc = orthopy.c1.jacobi.RecurrenceCoefficients("monic", alpha, beta, symbolic=True)
    _, a, b = np.array([rc[k] for k in range(n)]).T
    points, weights = scheme_from_rc(a, b, rc.int_1, mode=mode)
    return C1Scheme("Gauss-Jacobi", degree, weights, points)
