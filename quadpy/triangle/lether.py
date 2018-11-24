# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from ..line_segment import GaussLegendre


class Lether(object):
    """
    Frank G. Lether,
    Computation of double integrals over a triangle,
    Journal of Computational and Applied Mathematics,
    Volume 2, Issue 3, September 1976, Pages 219–224,
    <https://doi.org/10.1016/0771-050X(76)90008-5>.

    The results have later been reproduced in:

    H.T. Rathod, K.V. Nagaraja, B. Venkatesudu,
    Symmetric Gauss Legendre quadrature formulas for composite numerical integration
    over a triangular surface,
    Applied Mathematics and Computation 188 (2007) 865–876
    <https://doi.org/10.1016/j.amc.2006.10.041>.
    """

    def __init__(self, n):
        self.name = "Lether({})".format(n)

        gl = GaussLegendre(n)

        w = numpy.outer((1 + gl.points) * gl.weights, gl.weights) / 4
        x = numpy.outer(1 - gl.points, numpy.ones(n)) / 2
        y = numpy.outer(1 + gl.points, 1 - gl.points) / 4

        self.points = numpy.array([x.flatten(), y.flatten()]).T
        self.weights = w.flatten()

        self.bary = numpy.array(
            [self.points[:, 0], self.points[:, 1], 1 - numpy.sum(self.points, axis=1)]
        ).T

        self.degree = 2 * (n - 1)
        return
