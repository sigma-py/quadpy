# -*- coding: utf-8 -*-
#
from math import pi, sqrt

from ..helpers import untangle, pm


class AlbrechtCollatz(object):
    # Look up citation
    '''
    Albrecht-Collatz [1]
    '''
    def __init__(self):
        # Stroud claims degree 3, but really it's only order 1.
        self.degree = 1
        data = [
            (0.25, pm(2, sqrt(0.5))),
            ]
        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
