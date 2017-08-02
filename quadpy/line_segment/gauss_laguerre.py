# -*- coding: utf-8 -*-
#
import numpy


class GaussLaguerre(object):
    '''
    Gauss-Laguerre quadrature.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.laguerre.laggauss(n)
        return
