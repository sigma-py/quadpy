# -*- coding: utf-8 -*-
#
import numpy


class GaussHermite(object):
    '''
    Gauss-Hermite quadrature for integrals of the form

        int_{-inf}^{+inf} exp(-x^2) f(x) dx.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.hermite.hermgauss(n)
        return
