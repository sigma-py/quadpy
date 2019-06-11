# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

import warnings


class CohenGismalla(object):
    """
     A. M. Cohen, D. A. Gismalla,
     Some integration formulae for symmetric functions of two variables,
     International Journal of Computer Mathematics, 1986, 19:1, 57-68,
     <https://doi.org/10.1080/00207168608803504>.
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        # TODO improve precision
        warnings.warn("The Cohen-Gismalla schemes are only given in single-precision.")

        self.name = "CohenGismalla({})".format(index)

        if index == 1:
            # This scheme is of order 5 for symmetric integrands
            self.degree = 3
            u = 0.84623312
            v = 0.46607171
            self.points = numpy.array(
                [[0.0, 0.0], [+u, -v], [-u, +v], [-v, -u], [+v, +u]]
            )
            self.weights = numpy.array(
                [frac(8, 7), frac(5, 7), frac(5, 7), frac(5, 7), frac(5, 7)]
            )
        else:
            assert index == 2

            # ERR this scheme only has order 1
            # According to the article, it has order 7 for symmetric integrands.
            # Something is fishy...
            self.degree = 1
            r = 0.5878606
            s = 0.9353943
            u = 0.6105540
            v = 0.1109710
            A = 0.1856914
            B = 0.5951448
            C = 0.3584324
            self.points = numpy.array(
                [
                    [0.0, 0.0],
                    [+u, -v],
                    [-u, +v],
                    [-v, -u],
                    [+v, +u],
                    [+r, -s],
                    [-r, +s],
                    [-r, -s],
                    [+r, +s],
                ]
            )
            self.weights = numpy.array([A, B, B, B, B, C, C, C, C])
        return
