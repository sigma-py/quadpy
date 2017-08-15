# -*- coding: utf-8 -*-
#
from math import pi, sqrt

from ..helpers import untangle, fsd


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
        else:
            assert False

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
