# -*- coding: utf-8 -*-
#
import orthopy

from ..tools import scheme_from_rc


class GaussJacobi(object):
    """
    Gauss-Jacobi quadrature.
    """

    def __init__(self, n, alpha, beta, mode="numpy", decimal_places=None):

        self.degree = 2 * n - 1

        _, _, a, b = orthopy.line_segment.recurrence_coefficients.jacobi(
            n, alpha, beta, "monic", symbolic=True
        )
        self.points, self.weights = scheme_from_rc(
            a, b, mode=mode, decimal_places=decimal_places
        )
        return
