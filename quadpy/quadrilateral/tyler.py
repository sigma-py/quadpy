# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import zero, symm_s, symm_r0, concat


class Tyler(object):
    """
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://doi.org/10.4153/CJM-1953-044-1>.
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        pm = numpy.array([+1, -1])

        self.name = "Tyler({})".format(index)
        if index == 1:
            self.degree = 5
            self.points, self.weights = concat(
                zero(-frac(28, 45)),
                symm_s([frac(1, 36), 1]),
                symm_r0([frac(1, 45), 1], [frac(16, 45), frac(1, 2)]),
            )
        elif index == 2:
            self.degree = 7
            r = sqrt(frac(6, 7))
            s, t = sqrt((114 - pm * 3 * sqrt(583)) / 287)
            B1 = frac(49, 810)
            B2, B3 = (178981 + pm * 2769 * sqrt(583)) / 1888920
            self.points, self.weights = concat(
                symm_r0([B1, r]), symm_s([B2, s], [B3, t])
            )
        else:
            assert index == 3
            self.degree = 7
            self.points, self.weights = concat(
                zero(frac(449, 315)),
                symm_r0(
                    [frac(37, 1260), 1],
                    [frac(3, 28), frac(2, 3)],
                    [-frac(69, 140), frac(1, 3)],
                ),
                symm_s([frac(7, 540), 1], [frac(32, 135), frac(1, 2)]),
            )

        self.weights *= 4
        return
