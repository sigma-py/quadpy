# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s, _z

from ..helpers import untangle


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
            data = [
                (-28.0/45.0, _z()),
                (1.0/36.0, _symm_s(1.0)),
                (1.0/45.0, _symm_r_0(1.0)),
                (16.0/45.0, _symm_r_0(0.5)),
                ]
        elif index == 2:
            self.degree = 7
            r = numpy.sqrt(6.0 / 7.0)
            s = numpy.sqrt((114.0 - 3*numpy.sqrt(583.0)) / 287.0)
            t = numpy.sqrt((114.0 + 3*numpy.sqrt(583.0)) / 287.0)
            B1 = 49.0 / 810.0
            B2 = (178981.0 + 2769 * numpy.sqrt(583.0)) / 1888920.0
            B3 = (178981.0 - 2769 * numpy.sqrt(583.0)) / 1888920.0
            data = [
                (B1, _symm_r_0(r)),
                (B2, _symm_s(s)),
                (B3, _symm_s(t)),
                ]
        else:
            assert index == 3
            self.degree = 7
            r = 2.0/3.0
            s = 1.0/3.0
            t = 0.5
            data = [
                (449.0/315.0, _z()),
                (37.0/1260.0, _symm_r_0(1.0)),
                (3.0/28.0, _symm_r_0(r)),
                (-69.0/140.0, _symm_r_0(s)),
                (7.0/540.0, _symm_s(1.0)),
                (32.0/135.0, _symm_s(t)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 4.0
        return
