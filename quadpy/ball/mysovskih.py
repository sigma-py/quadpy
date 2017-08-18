# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle, fsd, pm


class Mysovskih(object):
    '''
    I.P. Mysovskih,
    On the construction of cubature formulas for the simplest regions,
    Z. Vychisl. Mat. i. Mat. Fiz. 4, 3-14, 1964.
    '''
    # pylint: disable=too-many-locals
    def __init__(self):
        self.degree = 7

        plus_minus = numpy.array([+1, -1])

        sqrt17770 = numpy.sqrt(17770.0)
        r, s = numpy.sqrt((1715.0 - plus_minus * 7 * sqrt17770) / 2817.0)
        t = numpy.sqrt(7.0/18.0)
        u = numpy.sqrt(7.0/27.0)

        B1, B2 = \
            (2965 * sqrt17770 + plus_minus * 227816) / 72030.0 / sqrt17770
        B3 = 324.0 / 12005.0
        B4 = 2187.0 / 96040.0

        data = [
            (B1, fsd(3, r, 1)),
            (B2, fsd(3, s, 1)),
            (B3, fsd(3, t, 2)),
            (B4, pm(3, u)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4.0/3.0 * numpy.pi
        return
