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

    def __init__(self, index, lmbda=None, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "Franke({})".format(index)

        if index == "1":
            assert lmbda is not None, "Franke(1) needs a lmbda"
            assert -frac(9, 5) <= lmbda <= frac(9, 4)

            self.degree = 5

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
        else:
            assert index == "2a"
            self.degree = 7

            a = sqrt(frac(15 + 2 * sqrt(30), 35))
            b = sqrt(frac(15 - 2 * sqrt(30), 35))

            self.points = numpy.array(
                [
                    [+0.105784012371275e1, +a],
                    [+0.105784012371275e1, -a],
                    [-0.105784012371275e1, +a],
                    [-0.105784012371275e1, -a],
                    [+0.774596669241483, +b],
                    [+0.774596669241483, -b],
                    [-0.774596669241483, +b],
                    [-0.774596669241483, -b],
                    [+0.469253522127911, +a],
                    [+0.469253522127911, -a],
                    [-0.469253522127911, +a],
                    [-0.469253522127911, -a],
                    [0, +b],
                    [0, -b],
                ]
            )
            self.weights = numpy.array(
                [
                    0.437841520872291e-1,
                    0.437841520872291e-1,
                    0.437841520872291e-1,
                    0.437841520872291e-1,
                    0.362302863812526,
                    0.362302863812526,
                    0.362302863812526,
                    0.362302863812526,
                    0.304070693050225,
                    0.304070693050225,
                    0.304070693050225,
                    0.304070693050225,
                    0.579684582100041,
                    0.579684582100041,
                ]
            )
        return
