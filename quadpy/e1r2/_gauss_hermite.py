import numpy

import orthopy

from ..tools import scheme_from_rc
from ._helpers import E1r2Scheme


def gauss_hermite(n, mode="numpy"):
    """
    Gauss-Hermite quadrature for integrals of the form

        int_{-inf}^{+inf} exp(-x^2) f(x) dx.
    """
    if mode == "numpy":
        points, weights = numpy.polynomial.hermite.hermgauss(n)
    else:
        _, _, alpha, beta = orthopy.e1r2.recurrence_coefficients(
            n, "monic", symbolic=True
        )
        # For some reason, the parameters have to be adapted here. TODO find out why
        beta[1:] /= 2
        points, weights = scheme_from_rc(alpha, beta, mode=mode)
    return E1r2Scheme(f"Gauss-Hermite ({n})", weights, points, 2 * n - 1)
