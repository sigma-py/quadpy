# -*- coding: utf-8 -*-
#
import math

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
                (1.0, fsd(2, math.sqrt(2.0/3.0), 1))
                ]
        elif index == '2-2':
            self.degree = 5
            alpha = math.sqrt(3.0/5.0)
            data = [
                (64.0/81.0, z(2)),
                (40.0/81.0, fsd(2, alpha, 1)),
                (25.0/81.0, pm(2, alpha)),
                ]
        else:
            assert index == '3-2'
            self.degree = 7
            alpha = math.sqrt(3.0/5.0)
            xi1 = math.sqrt(3.0/287.0 * (38.0 - math.sqrt(583.0)))
            xi2 = math.sqrt(3.0/287.0 * (38.0 + math.sqrt(583.0)))
            data = [
                (98.0/405.0, fsd(2, math.sqrt(6.0/7.0), 1)),
                (0.5205929166673945, pm(2, xi1)),
                (0.2374317746906302, pm(2, xi2)),
                ]

        self.points, self.weights = untangle(data)
        return
