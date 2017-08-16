# -*- coding: utf-8 -*-
#
from math import sqrt, pi
import numpy

from ..helpers import untangle, z, fsd, fs_array


class Mysovskih(object):
    '''
    I.P. Mysovskih,
    On the construction of cubature formulas for the simplest regions,
    Z. Vychisl. Mat. i. Mat. Fiz. 4, 3-14, 1964.
    '''
    def __init__(self, index, alpha=0.0):
        if index == 1:
            self.degree = 4
            i = numpy.arange(5)
            b = sqrt((alpha + 4.0)/(alpha + 6.0))

            r = b * numpy.cos(0.4*i*pi)
            s = b * numpy.sin(0.4*i*pi)
            x = numpy.array([r, s]).T

            B0 = 4.0 / (alpha + 4)**2
            B1 = (alpha + 2.0)*(alpha + 6.0) / 5.0 / (alpha + 4.0)**2

            data = [
                (B0, z(2)),
                (B1, x),
                ]
        else:
            assert index == 2
            self.degree = 11

            sqrt10 = sqrt(10.0)
            sqrt601 = sqrt(601.0)

            B1 = (857.0*sqrt601 + 12707.0) / 20736.0 / sqrt601
            B2 = 125.0 / 3456.0
            B3 = (857.0*sqrt601 - 12707.0) / 20736.0 / sqrt601
            B4 = (340.0 + 25*sqrt10) / 10368.0
            B5 = (340.0 - 25*sqrt10) / 10368.0

            r1 = sqrt((31.0 - sqrt601) / 60.0)
            r2 = sqrt(3.0/5.0)
            r3 = sqrt((31.0 + sqrt601) / 60.0)
            r4 = sqrt((10.0 - sqrt10) / 20.0)
            r5 = sqrt((10.0 + sqrt10) / 20.0)

            s4 = sqrt((10.0 - sqrt10) / 60.0)
            s5 = sqrt((10.0 + sqrt10) / 60.0)

            data = [
                (B1, fsd(2, r1, 1)),
                (B2, fsd(2, r2, 1)),
                (B3, fsd(2, r3, 1)),
                (B4, fs_array([r4, s4])),
                (B5, fs_array([r5, s5])),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= pi
        return
