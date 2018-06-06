# -*- coding: utf-8 -*-
#
import numpy
import orthopy

from ..tools import scheme_from_rc


class GaussHermite(object):
    """
    Gauss-Hermite quadrature for integrals of the form

        int_{-inf}^{+inf} exp(-x^2) f(x) dx.
    """

    def __init__(self, n, mode="numpy", decimal_places=None):
        self.degree = 2 * n - 1

        if mode == "numpy":
            self.points, self.weights = numpy.polynomial.hermite.hermgauss(n)
        else:
            _, _, alpha, beta = orthopy.e1r2.recurrence_coefficients(
                n, "monic", symbolic=True
            )

            # For some reason, the parameters have to be adapted here.
            # TODO find out why
            beta[1:] /= 2

            self.points, self.weights = scheme_from_rc(
                alpha, beta, mode=mode, decimal_places=decimal_places
            )
        return
