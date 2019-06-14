# -*- coding: utf-8 -*-
#
import orthopy

from .helpers import E1rScheme
from ..tools import scheme_from_rc


def GaussLaguerre(n, alpha=0, mode="numpy", decimal_places=None):
    """
    Gauss-Laguerre quadrature for integrals of the form

        int_0^{+inf} exp(-x) f(x) dx.
    """
    _, _, a, b = orthopy.e1r.recurrence_coefficients(n, alpha, "monic", symbolic=True)
    points, weights = scheme_from_rc(a, b, mode=mode, decimal_places=decimal_places)
    return E1rScheme("Gauss-Laguerre", 2 * n - 1, weights, points)
