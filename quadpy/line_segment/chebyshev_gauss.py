# -*- coding: utf-8 -*-
#
import numpy


class ChebyshevGauss1(object):
    '''
    Chebyshev-Gauss quadrature for \\int_{-1}^1 f(x) / sqrt(1+x^2) dx.
    '''
    def __init__(self, n):
        self.degree = n if n % 2 == 1 else n+1
        self.points = numpy.cos(
                (2*numpy.arange(1, n+1) - 1.0) / (2*n) * numpy.pi
                )
        self.weights = numpy.full(n, numpy.pi / n)
        return


class ChebyshevGauss2(object):
    '''
    Chebyshev-Gauss quadrature for \\int_{-1}^1 f(x) * sqrt(1+x^2) dx.
    '''
    def __init__(self, n):
        self.degree = n if n % 2 == 1 else n+1
        self.points = numpy.cos(
                numpy.pi * numpy.arange(1, n+1) / (n+1)
                )
        self.weights = numpy.pi / (n+1) \
            * (numpy.sin(numpy.pi * numpy.arange(1, n+1) / (n+1)))**2
        return
