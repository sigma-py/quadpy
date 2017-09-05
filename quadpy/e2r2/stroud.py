# -*- coding: utf-8 -*-
#
from __future__ import division

import warnings

import numpy

from .rabinowitz_richter import RabinowitzRichter
from .stroud_secrest import StroudSecrest

from ..helpers import untangle, fsd, pm


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    def __init__(self, index):
        self.name = 'Stroud_E2r2({})'.format(index)
        if index == '4-1':
            self.degree = 4

            i = numpy.arange(5)
            pts = numpy.sqrt(2.0) * numpy.array([
                numpy.cos(2*i*numpy.pi / 5),
                numpy.sin(2*i*numpy.pi / 5)
                ]).T
            data = [
                (0.5, numpy.array([[0.0, 0.0]])),
                (0.1, pts),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= numpy.pi
        elif index == '5-1':
            self.set_data(StroudSecrest('V'))
        elif index == '5-2':
            # Cartesian product Gauss formula
            self.degree = 5

            r = numpy.sqrt(3/2)
            data = [
                (4/9, numpy.array([[0.0, 0.0]])),
                (1/9, fsd(2, (r, 1))),
                (1/36, pm(2, r)),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= numpy.pi
        elif index == '7-1':
            self.set_data(StroudSecrest('VI'))
        elif index == '7-2':
            # Cartesian product Gauss formula
            warnings.warn(
                'Stroud\'s Gauss product formula has degree 1, not 7.'
                )
            self.degree = 1

            p_m = numpy.array([+1, -1])
            sqrt6 = numpy.sqrt(6)
            r, s = numpy.sqrt((3 + p_m * sqrt6) / 2)
            A, B = (5 - p_m * 2 * sqrt6) / 48
            C = 1/48

            data = [
                (A, fsd(2, (r, 1))),
                (B, fsd(2, (s, 1))),
                (C, fsd(2, (r, 1), (s, 1))),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= numpy.pi
        elif index == '9-1':
            self.set_data(RabinowitzRichter(1))
        elif index == '11-1':
            self.set_data(RabinowitzRichter(2))
        elif index == '11-2':
            self.set_data(RabinowitzRichter(3))
        elif index == '13-1':
            self.set_data(RabinowitzRichter(4))
        else:
            assert index == '15-1'
            self.set_data(RabinowitzRichter(5))
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
