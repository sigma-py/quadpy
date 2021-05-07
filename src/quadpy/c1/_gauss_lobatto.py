import numpy as np
import orthopy

from ..tools import scheme_from_rc
from ._helpers import C1Scheme


def gauss_lobatto(n, a=0.0, b=0.0):
    assert n >= 2
    degree = 2 * n - 3
    rc = orthopy.c1.jacobi.RecurrenceCoefficients("monic", a, b, symbolic=False)
    _, alpha, beta = np.array([rc[k] for k in range(n)]).T
    beta[0] = rc.int_1
    points, weights = _lobatto(alpha, beta, rc.int_1, -1.0, 1.0)
    return C1Scheme("Gauss-Lobatto", degree, weights, points)


def _lobatto(alpha, beta, int_1, xl1, xl2):
    """Compute the Lobatto nodes and weights with the preassigned node xl1, xl2.
    Based on the section 7 of the paper

        Some modified matrix eigenvalue problems,
        Gene Golub,
        SIAM Review Vol 15, No. 2, April 1973, pp.318--334,

    and

        http://www.scientificpython.net/pyblog/radau-quadrature
    """
    from scipy.linalg import solve, solve_banded

    n = len(alpha) - 1
    en = np.zeros(n)
    en[-1] = 1
    A1 = np.vstack((np.sqrt(beta), alpha - xl1))
    J1 = np.vstack((A1[:, 0:-1], A1[0, 1:]))
    A2 = np.vstack((np.sqrt(beta), alpha - xl2))
    J2 = np.vstack((A2[:, 0:-1], A2[0, 1:]))
    g1 = solve_banded((1, 1), J1, en)
    g2 = solve_banded((1, 1), J2, en)
    C = np.array(((1, -g1[-1]), (1, -g2[-1])))
    xl = np.array((xl1, xl2))
    ab = solve(C, xl)

    alphal = alpha.copy()
    alphal[-1] = ab[0]
    betal = beta.copy()
    betal[-1] = ab[1]
    x, w = scheme_from_rc(alphal, betal, int_1, mode="numpy")
    return x, w
