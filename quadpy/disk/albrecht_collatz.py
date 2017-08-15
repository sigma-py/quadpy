# -*- coding: utf-8 -*-
#
from math import pi

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
        self.degree = 3
        data = [
            # ERR Wrongly stated in Stroud as sqrt(0.5) instead of 0.5
            (0.25, pm(2, 0.5)),
            ]
        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
