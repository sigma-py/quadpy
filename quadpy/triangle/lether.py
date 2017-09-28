# -*- coding: utf-8 -*-
#
from __future__ import division

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

        w = numpy.outer((1 + a) * A, A)
        x = numpy.outer(1-a, numpy.ones(n)) / 2
        y = numpy.outer(1+a, 1-a) / 4

        self.weights = w.reshape(-1) / 4
        self.points = numpy.stack([x.reshape(-1), y.reshape(-1)]).T

        self.bary = numpy.array([
            self.points[:, 0],
            self.points[:, 1],
            1 - numpy.sum(self.points, axis=1)
            ]).T
        return
