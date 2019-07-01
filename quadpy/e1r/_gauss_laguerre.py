import orthopy

from ..tools import scheme_from_rc
from ._helpers import E1rScheme


def gauss_laguerre(n, alpha=0, mode="numpy"):
    """
    Gauss-Laguerre quadrature for integrals of the form

        int_0^{+inf} exp(-x) f(x) dx.
    """
    _, _, a, b = orthopy.e1r.recurrence_coefficients(n, alpha, "monic", symbolic=True)
    points, weights = scheme_from_rc(a, b, mode=mode)

    if alpha != 0:
        name = "Generalized Gauss-Laguerre (n={}, alpha={})".format(n, alpha)
    else:
        name = "Gauss-Laguerre ({})".format(n, alpha)

    return E1rScheme(name, weights, points, 2 * n - 1)
