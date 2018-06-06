# -*- coding: utf-8 -*-
#
import numpy
import orthopy

from ..tools import scheme_from_rc


class GaussLegendre(object):
    """
    Gauss-Legendre quadrature.
    """

    def __init__(self, n, mode="numpy", decimal_places=None):

        self.degree = 2 * n - 1

        if mode == "numpy":
            self.points, self.weights = numpy.polynomial.legendre.leggauss(n)
        else:
            _, _, alpha, beta = orthopy.line_segment.recurrence_coefficients.legendre(
                n, "monic", symbolic=True
            )
            self.points, self.weights = scheme_from_rc(
                alpha, beta, mode=mode, decimal_places=decimal_places
            )
        return
