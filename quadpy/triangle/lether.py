# -*- coding: utf-8 -*-
#
import numpy


class Lether(object):
    '''
    Frank G. Lether,
    Computation of double integrals over a triangle,
    Journal of Computational and Applied Mathematics,
    Volume 2, Issue 3, September 1976, Pages 219–224,
    <https://doi.org/10.1016/0771-050X(76)90008-5>.
    '''
    def __init__(self, n):
        self.degree = 2*n - 2

        a, A = numpy.polynomial.legendre.leggauss(n)

        w = numpy.outer((1.0 + a) * A, A)
        x = 0.50 * numpy.outer(1.0-a, numpy.ones(n))
        y = 0.25 * numpy.outer(1+a, 1-a)

        self.weights = 0.25 * w.reshape(-1)
        self.points = numpy.stack([x.reshape(-1), y.reshape(-1)]).T

        self.bary = numpy.array([
            self.points[:, 0],
            self.points[:, 1],
            1.0 - numpy.sum(self.points, axis=1)
            ]).T

        return
