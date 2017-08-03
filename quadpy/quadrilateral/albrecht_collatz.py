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
            self.weights = numpy.concatenate([
                numpy.full(1, 5.0/12.0),
                numpy.full(4, 0.125),
                numpy.full(4, 1.0/48.0),
                ])
            self.points = numpy.concatenate([
                _z(),
                _symm_r_0(1.0),
                _symm_s(1.0)
                ])
        elif index == 2:
            self.degree = 5
            self.weights = numpy.concatenate([
                numpy.full(4, 5.0/36.0),
                numpy.full(2, 5.0/63.0),
                numpy.full(1, 2.0/7.0),
                ])
            r = numpy.sqrt(3.0 / 5.0)
            s = numpy.sqrt(1.0 / 3.0)
            t = numpy.sqrt(14.0 / 15.0)
            self.points = numpy.concatenate([
                _pm2(r, s),
                _pm(0.0, t),
                _z(),
                ])
        elif index == 3:
            self.degree = 5
            self.weights = numpy.concatenate([
                numpy.full(1, 2.0/7.0),
                numpy.full(2, 25.0/168.0),
                numpy.full(2, 5.0/48.0),
                numpy.full(2, 5.0/48.0),
                ])
            r = numpy.sqrt(7.0 / 15.0)
            s = numpy.sqrt((7.0 + numpy.sqrt(24)) / 15.0)
            t = numpy.sqrt((7.0 - numpy.sqrt(24)) / 15.0)
            self.points = numpy.concatenate([
                _z(),
                _pm(r, r),
                _pm(+s, -t),
                _pm(+t, -s),
                ])
        else:
            assert index == 4
            self.degree = 5
            self.weights = numpy.concatenate([
                numpy.full(1, 2.0/45.0),
                numpy.full(4, 2.0/45.0),
                numpy.full(4, 1.0/60.0),
                numpy.full(4, 8.0/45.0),
                ])
            self.points = numpy.concatenate([
                _z(),
                _symm_r_0(1.0),
                _symm_s(1.0),
                _symm_s(0.5),
                ])

        self.weights *= 4.0
        return
