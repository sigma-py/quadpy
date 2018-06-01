# -*- coding: utf-8 -*-
#
import numpy


class Lether(object):
    """
    Frank G. Lether,
    A Generalized Product Rule for the Circle,
    SIAM Journal on Numerical Analysis,
    Vol. 8, No. 2 (Jun., 1971), pp. 249-253,
    <http://www.jstor.org/stable/2949473>.
    """

    def __init__(self, n):
        p, w = numpy.polynomial.legendre.leggauss(n)

        mu = numpy.arange(1, n + 1)
        self.points = numpy.column_stack(
            [
                numpy.tile(numpy.cos(mu * numpy.pi / (n + 1)), n),
                numpy.outer(p, numpy.sin(mu * numpy.pi / (n + 1))).flatten(),
            ]
        )

        self.weights = (
            numpy.pi
            / (n + 1)
            * numpy.outer(w, numpy.sin(mu * numpy.pi / (n + 1)) ** 2).flatten()
        )

        self.degree = 2 * n - 1
        return
