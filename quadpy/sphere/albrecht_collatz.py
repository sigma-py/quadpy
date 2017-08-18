# -*- coding: utf-8 -*-
#
import math
import numpy

from ..helpers import untangle, pm_array0, fsd, pm


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://dx.doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self, index):
        if index == 1:
            self.degree = 5
            plus_minus = numpy.array([+1, -1])
            r, s = numpy.sqrt((5 + plus_minus * math.sqrt(5.0)) / 10.0)
            data = [
                (1.0/12.0, pm_array0(3, [r, s], [0, 1])),
                (1.0/12.0, pm_array0(3, [r, s], [1, 2])),
                (1.0/12.0, pm_array0(3, [r, s], [2, 0])),
                ]
        elif index == 2:
            self.degree = 5
            r = 1.0
            s = math.sqrt(1.0/3.0)
            data = [
                (8.0/120.0, fsd(3, r, 1)),
                (9.0/120.0, pm(3, s)),
                ]
        elif index == 3:
            self.degree = 5
            r = 1.0
            s = math.sqrt(1.0/2.0)
            data = [
                (1.0/30.0, fsd(3, r, 1)),
                (2.0/30.0, fsd(3, s, 2)),
                ]
        elif index == 4:
            self.degree = 5
            plus_minus = numpy.array([+1, -1])
            r, s = numpy.sqrt((3 + plus_minus * math.sqrt(5.0)) / 6.0)
            t = math.sqrt(1.0/3.0)
            data = [
                (1.0/20.0, pm_array0(3, [r, s], [0, 1])),
                (1.0/20.0, pm_array0(3, [r, s], [1, 2])),
                (1.0/20.0, pm_array0(3, [r, s], [2, 0])),
                (1.0/20.0, pm(3, t)),
                ]
        else:
            assert False

        self.points, self.weights = untangle(data)
        return
