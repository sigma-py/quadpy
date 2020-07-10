import numpy
import orthopy

from ..tools import scheme_from_rc
from ._helpers import E1rScheme


def gauss_laguerre(n, alpha=0, mode="numpy"):
    """
    Gauss-Laguerre quadrature for integrals of the form

        int_0^{+inf} exp(-x) f(x) dx.
    """
    symbolic = (mode != "numpy")
    rc = orthopy.e1r.RecurrenceCoefficients("monic", alpha, symbolic)
    _, a, b = numpy.array([rc[k] for k in range(n)]).T
    b[0] = rc.int_1

    points, weights = scheme_from_rc(a, b, mode=mode)

    if alpha == 0:
        name = f"Gauss-Laguerre ({n})"
    else:
        name = f"Generalized Gauss-Laguerre (n={n}, alpha={alpha})"

    return E1rScheme(name, weights, points, 2 * n - 1)
