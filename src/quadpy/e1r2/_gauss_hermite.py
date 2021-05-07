import numpy as np
import orthopy

from ..tools import scheme_from_rc
from ._helpers import E1r2Scheme


def gauss_hermite(n, mode="numpy"):
    """
    Gauss-Hermite quadrature for integrals of the form

        int_{-inf}^{+inf} exp(-x^2) f(x) dx.
    """
    if mode == "numpy":
        points, weights = np.polynomial.hermite.hermgauss(n)
    else:
        rc = orthopy.e1r2.RecurrenceCoefficients("physicists", "monic", symbolic=True)
        _, alpha, beta = np.array([rc[k] for k in range(n)]).T
        points, weights = scheme_from_rc(alpha, beta, rc.int_1, mode=mode)
    return E1r2Scheme(f"Gauss-Hermite ({n})", weights, points, 2 * n - 1)
