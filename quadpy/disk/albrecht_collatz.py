# -*- coding: utf-8 -*-
#
from math import pi, sqrt
import warnings

from ..helpers import untangle, pm


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://dx.doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self):
        self.name = 'AlbrechtCollatz'
        warnings.warn('Formula {} only has degree 1, not 3!'.format(self.name))
        self.degree = 1
        data = [
            (0.25, pm(2, sqrt(0.5))),
            ]
        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
