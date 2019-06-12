# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy


class Franke(object):
    """
    Richard Franke,
    Obtaining cubatures for rectangles and other planar regions by using orthogonal
    polynomials,
    Math. Comp. 25 (1971), 803-817,
    <https://doi.org/10.1090/S0025-5718-1971-0300440-5>.
    """

    def __init__(self, lmbda, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 5

        assert -frac(9, 5) <= lmbda <= frac(9, 4)

        self.name = "Franke"

        a = sqrt(frac(9 + 5 * lmbda, 15))
        b = sqrt(frac(9 - 4 * lmbda, 15))
        c = sqrt(frac(3, 5))

        self.points = numpy.array(
            [
                [0, 0],
                [+a, 0],
                [-a, 0],
                [+b, +c],
                [+b, -c],
                [-b, +c],
                [-b, -c],
                [0, +c],
                [0, -c],
            ]
        )
        self.weights = numpy.array(
            [
                frac(16 * (4 + 5 * lmbda), 9 * (9 + 5 * lmbda)),
                frac(40, 9 * (9 + 5 * lmbda)),
                frac(40, 9 * (9 + 5 * lmbda)),
                frac(25, 9 * (9 - 4 * lmbda)),
                frac(25, 9 * (9 - 4 * lmbda)),
                frac(25, 9 * (9 - 4 * lmbda)),
                frac(25, 9 * (9 - 4 * lmbda)),
                frac(40 * (1 - lmbda), 9 * (9 - 4 * lmbda)),
                frac(40 * (1 - lmbda), 9 * (9 - 4 * lmbda)),
            ]
        )
        return
