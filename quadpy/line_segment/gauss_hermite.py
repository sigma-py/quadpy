# -*- coding: utf-8 -*-
#
import numpy


class GaussHermite(object):
    '''
    Gauss-Hermite quadrature.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.hermite.hermgauss(n)
        return
