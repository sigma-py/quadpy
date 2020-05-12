import orthopy

from ..tools import scheme_from_rc
from ._helpers import C1Scheme


def gauss_jacobi(n, alpha, beta, mode="numpy"):
    degree = 2 * n - 1

    _, _, a, b = orthopy.line_segment.recurrence_coefficients.jacobi(
        n, alpha, beta, "monic", symbolic=True
    )
    points, weights = scheme_from_rc(a, b, mode=mode)
    return C1Scheme("Gauss-Jacobi", degree, weights, points)
