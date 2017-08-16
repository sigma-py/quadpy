# -*- coding: utf-8 -*-
#
import math

from .helpers import _s4, _s31
from ..helpers import untangle


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, degree):
        self.degree = degree
        if degree == 2:
            data = [
                (1.0/4.0, _s31((5.0 - math.sqrt(5.0))/20.0)),
                ]
        else:
            assert degree == 3
            data = [
                (-4.0/5.0, _s4()),
                (+9.0/20.0, _s31(1.0/6.0)),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
