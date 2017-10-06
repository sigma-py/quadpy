# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from ..helpers import untangle, fsd, z, pm


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, index):

        if index == '1-2':
            self.degree = 3
            data = [
                (1, fsd(2, (sqrt(fr(2, 3)), 1)))
                ]
        elif index == '2-2':
            self.degree = 5
            alpha = sqrt(fr(3, 5))
            data = [
                (fr(64, 81), z(2)),
                (fr(40, 81), fsd(2, (alpha, 1))),
                (fr(25, 81), pm(2, alpha)),
                ]
        else:
            assert index == '3-2'
            self.degree = 7
            alpha = sqrt(fr(3, 5))
            xi1, xi2 = [
                sqrt(fr(3, 287) * (38 - i*sqrt(583)))
                for i in [+1, -1]
                ]
            data = [
                (fr(98, 405), fsd(2, (sqrt(fr(6, 7)), 1))),
                (0.5205929166673945, pm(2, xi1)),
                (0.2374317746906302, pm(2, xi2)),
                ]

        self.points, self.weights = untangle(data)
        return
