# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s, _z, _pm, _pm2


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://dx.doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self, index):
        self.name = 'AlbrechtCollatz({})'.format(index)
        if index == 1:
            self.degree = 3
            data = [
                (_z(), 5.0/12.0),
                (_symm_r_0(1.0), 0.125),
                (_symm_s(1.0), 1.0/48.0)
                ]
        elif index == 2:
            self.degree = 5
            r = numpy.sqrt(3.0 / 5.0)
            s = numpy.sqrt(1.0 / 3.0)
            t = numpy.sqrt(14.0 / 15.0)
            data = [
                (_pm2(r, s), 5.0/36.0),
                (_pm(0.0, t), 5.0/63.0),
                (_z(), 2.0/7.0)
                ]
        elif index == 3:
            self.degree = 5
            r = numpy.sqrt(7.0 / 15.0)
            s = numpy.sqrt((7.0 + numpy.sqrt(24)) / 15.0)
            t = numpy.sqrt((7.0 - numpy.sqrt(24)) / 15.0)
            data = [
                (_z(), 2.0/7.0),
                (_pm(r, r), 25.0/168.0),
                (_pm(+s, -t), 5.0/48.0),
                (_pm(+t, -s), 5.0/48.0),
                ]
        else:
            assert index == 4
            self.degree = 5
            data = [
                (_z(), 2.0/45.0),
                (_symm_r_0(1.0), 2.0/45.0),
                (_symm_s(1.0), 1.0/60.0),
                (_symm_s(0.5), 8.0/45.0),
                ]

        points, weights = zip(*data)
        self.points = numpy.concatenate(points)
        self.weights = numpy.repeat(
                weights, [len(group) for group in points]
                )
        self.weights *= 4.0
        return
