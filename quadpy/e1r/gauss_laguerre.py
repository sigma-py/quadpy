# -*- coding: utf-8 -*-
#
import orthopy

from ..tools import scheme_from_rc


class GaussLaguerre(object):
    '''
    Gauss-Laguerre quadrature for integrals of the form

        int_0^{+inf} exp(-x) f(x) dx.
    '''
    def __init__(self, n, alpha=0, mode='numpy', decimal_places=None):
        self.degree = 2*n - 1

        _, _, a, b = orthopy.e1r.recurrence_coefficients(
                n, alpha, 'monic', symbolic=True
                )
        self.points, self.weights = scheme_from_rc(
                a, b, mode=mode, decimal_places=decimal_places
                )
        return
