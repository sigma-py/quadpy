# -*- coding: utf-8 -*-
#
from math import pi, sqrt

from ..helpers import untangle, fsd, pm, z


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, index):
        if index == '11-2':
            self.degree = 3
            # ERR Wrongly stated in Stroud with 0.5 instead of sqrt(0.5)
            data = [
                (0.25, fsd(2, sqrt(0.5), 1)),
                ]
        elif index == '12-2':
            self.degree = 5
            data = [
                (1.0/6.0, z(2)),
                (1.0/6.0, fsd(2, sqrt(0.5), 1)),
                (1.0/24.0, pm(2, sqrt(0.5))),
                ]
        else:
            assert index == '13-2'
            self.degree = 7
            sqrt29 = sqrt(29.0)
            b1 = (551.0 + 41.0 * sqrt29) / 6264.0
            b2 = (551.0 - 41.0 * sqrt29) / 6264.0
            xi1 = sqrt(3.0 / 2 / (9 + sqrt29))
            xi2 = sqrt(3.0 / 2 / (9 - sqrt29))
            data = [
                (2.0/27.0, fsd(2, sqrt(0.75), 1)),
                (b1, pm(2, xi1)),
                (b2, pm(2, xi2)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
