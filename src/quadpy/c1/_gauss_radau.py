import numpy as np
import orthopy

from ..tools import scheme_from_rc
from ._helpers import C1Scheme


def gauss_radau(n, a=0.0, b=0.0):
    assert n >= 2
    degree = 2 * n - 1

    rc = orthopy.c1.jacobi.RecurrenceCoefficients("monic", a, b, symbolic=False)
    _, alpha, beta = np.array([rc[k] for k in range(n)]).T
    points, weights = _radau(alpha, beta, rc.int_1, -1.0)
    return C1Scheme("Gauss-Radau", degree, weights, points)


def _radau(alpha, beta, int_1, xr):
    """From <http://www.scientificpython.net/pyblog/radau-quadrature>:
    Compute the Radau nodes and weights with the preassigned node xr.

    Based on the section 7 of the paper

        Some modified matrix eigenvalue problems,
        Gene Golub,
        SIAM Review Vol 15, No. 2, April 1973, pp.318--334.
    """
    from scipy.linalg import solve_banded

    beta[0] = int_1

    n = len(alpha) - 1
    f = np.zeros(n)
    f[-1] = beta[-1]
    A = np.vstack((np.sqrt(beta), alpha - xr))
    J = np.vstack((A[:, 0:-1], A[0, 1:]))
    delta = solve_banded((1, 1), J, f)
    alphar = alpha.copy()
    alphar[-1] = xr + delta[-1]
    x, w = scheme_from_rc(alphar, beta, int_1, mode="numpy")
    return x, w
