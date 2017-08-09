# -*- coding: utf-8 -*-
#
from math import pi

from ..helpers import untangle, fsd


class HammerStroud(object):
    # Look up citation
    '''
    Hammer-Stroud [2]
    '''
    def __init__(self):
        # Stroud claims degree 3, but really it's only order 1.
        self.degree = 1
        data = [
            (0.25, fsd(2, 0.5, 1)),
            ]
        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
