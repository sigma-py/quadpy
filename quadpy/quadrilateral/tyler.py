# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s, _z


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://dx.doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self, index):
        self.name = 'Tyler({})'.format(index)
        if index == 1:
            self.degree = 5
            self.weights = numpy.concatenate([
                numpy.full(1, -28.0/45.0),
                numpy.full(4, 1.0/36.0),
                numpy.full(4, 1.0/45.0),
                numpy.full(4, 16.0/45.0),
                ])
            self.points = numpy.concatenate([
                _z(),
                _symm_s(1.0),
                _symm_r_0(1.0),
                _symm_r_0(0.5),
                ])
        elif index == 2:
            self.degree = 7
            r = numpy.sqrt(6.0 / 7.0)
            s = numpy.sqrt((114.0 - 3*numpy.sqrt(583.0)) / 287.0)
            t = numpy.sqrt((114.0 + 3*numpy.sqrt(583.0)) / 287.0)
            B1 = 49.0 / 810.0
            B2 = (178981.0 + 2769 * numpy.sqrt(583.0)) / 1888920.0
            B3 = (178981.0 - 2769 * numpy.sqrt(583.0)) / 1888920.0
            self.weights = numpy.concatenate([
                numpy.full(4, B1),
                numpy.full(4, B2),
                numpy.full(4, B3),
                ])
            self.points = numpy.concatenate([
                _symm_r_0(r),
                _symm_s(s),
                _symm_s(t),
                ])
        else:
            assert index == 3
            self.degree = 7
            r = 2.0/3.0
            s = 1.0/3.0
            t = 0.5
            self.weights = numpy.concatenate([
                numpy.full(1, 449.0/315.0),
                numpy.full(4, 37.0/1260.0),
                numpy.full(4, 3.0/28.0),
                numpy.full(4, -69.0/140.0),
                numpy.full(4, 7.0/540.0),
                numpy.full(4, 32.0/135.0),
                ])
            self.points = numpy.concatenate([
                _z(),
                _symm_r_0(1.0),
                _symm_r_0(r),
                _symm_r_0(s),
                _symm_s(1.0),
                _symm_s(t),
                ])

        self.weights *= 4.0
        return
