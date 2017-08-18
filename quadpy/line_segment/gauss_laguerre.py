# -*- coding: utf-8 -*-
#
import numpy


class GaussLaguerre(object):
    '''
    Gauss-Laguerre quadrature for integrals of the form

        int_{-inf}^{+inf} exp(-x) f(x) dx.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.laguerre.laggauss(n)
        return
