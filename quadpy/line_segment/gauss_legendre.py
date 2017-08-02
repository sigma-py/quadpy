# -*- coding: utf-8 -*-
#
import numpy


class GaussLegendre(object):
    '''
    Gauss-Legendre quadrature.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.legendre.leggauss(n)
        return
