# -*- coding: utf-8 -*-
#
import numpy
import sympy

from ..helpers import untangle, fsd, pm


class Mysovskih(object):
    '''
    I.P. Mysovskih,
    On the construction of cubature formulas for the simplest regions,
    Z. Vychisl. Mat. i. Mat. Fiz. 4, 3-14, 1964.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        pi = sympy.pi if symbolic else numpy.pi
        pm_ = numpy.array([+1, -1])

        self.degree = 7

        sqrt17770 = sqrt(17770)
        r, s = sqrt((1715 - pm_ * 7 * sqrt17770) / 2817)
        t = sqrt(frac(7, 18))
        u = sqrt(frac(7, 27))

        B1, B2 = (2965 * sqrt17770 + pm_ * 227816) / 72030 / sqrt17770
        B3 = frac(324, 12005)
        B4 = frac(2187, 96040)

        data = [
            (B1, fsd(3, (r, 1))),
            (B2, fsd(3, (s, 1))),
            (B3, fsd(3, (t, 2))),
            (B4, pm(3, u)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= frac(4, 3) * pi
        return
