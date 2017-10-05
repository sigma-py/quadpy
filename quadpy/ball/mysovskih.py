# -*- coding: utf-8 -*-
#
from sympy import sqrt, pi, Rational as fr

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

        sqrt17770 = sqrt(17770)
        r, s = [
            sqrt((1715 - plus_minus * 7 * sqrt17770) / 2817)
            for plus_minus in [+1, -1]
            ]
        t = sqrt(fr(7, 18))
        u = sqrt(fr(7, 27))

        B1, B2 = [
            (2965 * sqrt17770 + plus_minus * 227816) / 72030 / sqrt17770
            for plus_minus in [+1, -1]
            ]
        B3 = fr(324, 12005)
        B4 = fr(2187, 96040)

        data = [
            (B1, fsd(3, (r, 1))),
            (B2, fsd(3, (s, 1))),
            (B3, fsd(3, (t, 2))),
            (B4, pm(3, u)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= fr(4, 3) * pi
        return
