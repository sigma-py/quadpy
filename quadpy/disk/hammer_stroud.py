# -*- coding: utf-8 -*-
#
from math import pi
import warnings

from ..helpers import untangle, fsd


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self):
        self.name = 'HammerStroud'
        warnings.warn('Formula {} only has degree 1, not 3!'.format(self.name))
        self.degree = 1
        data = [
            (0.25, fsd(2, 0.5, 1)),
            ]
        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
