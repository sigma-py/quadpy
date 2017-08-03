# -*- coding: utf-8 -*-
#
import numpy

from .miller import Miller
from .helpers import _symm_r_0, _symm_s, _symm_s_t

from .. import line_segment


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.

    <https://people.sc.fsu.edu/~jburkardt/m_src/stroud/square_unit_set.m>
    <http://nines.cs.kuleuven.be/ecf/mtables.html>
    '''
    def __init__(self, index):
        if index == 'C2 1-1':
            # Product trapezoidal
            self.degree = 1
            self.weights = numpy.full(4, 0.25)
            self.points = _symm_s(1.0)
        elif index == 'C2 1-2':
            m = Miller()
            self.degree = m.degree
            self.weights = m.weights
            self.points = m.points
        elif index == 'C2 3-1':
            # Product Gauss
            self.degree = 3
            self.weights = numpy.full(4, 0.25)
            self.points = _symm_s(1.0/numpy.sqrt(3.0))
        elif index == 2:
            self.weights = 4 * [1.0]
            self.points = _symm_r_0(numpy.sqrt(2.0/3.0))
            self.degree = 3
        elif index == 3:
            self.weights = (
                4 * [2.0/3.0] +
                [4.0/3.0]
                )
            self.points = (
                _symm_r_0(1.0) +
                [[0.0, 0.0]]
                )
            self.degree = 3
        elif index == 4:
            self.weights = (
                4 * [1.0/3.0] +
                [8.0/3.0]
                )
            self.points = (
                _symm_s(1.0) +
                [[0.0, 0.0]]
                )
            self.degree = 3
        elif index == 5:
            self.weights = (
                4 * [40.0 / 49.0] +
                4 * [9.0 / 49.0]
                )
            self.points = (
                _symm_r_0(numpy.sqrt(7.0/15.0)) +
                _symm_s(numpy.sqrt(7.0/9.0))
                )
            self.degree = 5
        elif index == 6:
            self.weights = (
                4 * [1.0 / 9.0] +
                1 * [-8.0 / 9.0] +
                4 * [10.0 / 9.0]
                )
            self.points = (
                _symm_s(1.0) +
                [[0.0, 0.0]] +
                _symm_r_0(numpy.sqrt(2.0 / 5.0))
                )
            self.degree = 5
        elif index == 7:
            self.weights = (
                [8.0 / 7.0] +
                8 * [5.0 / 14.0]
                )
            self.points = (
                [[0.0, 0.0]] +
                _symm_s_t(
                    0.846233119448533574334773553209421,
                    0.466071712166418914664132375714293
                    )
                )
            self.degree = 5
        elif index == 8:
            # Stroud number C2:7-1.
            r = numpy.sqrt(6.0 / 7.0)
            c = 3.0 * numpy.sqrt(583.0)
            s = numpy.sqrt((114.0 - c) / 287.0)
            t = numpy.sqrt((114.0 + c) / 287.0)
            w1 = 4.0 * 49.0 / 810.0
            w2 = 4.0 * (178981.0 + 923.0 * c) / 1888920.0
            w3 = 4.0 * (178981.0 - 923.0 * c) / 1888920.0
            #
            self.weights = (
                4 * [w1] +
                4 * [w2] +
                4 * [w3]
                )
            self.points = (
                _symm_r_0(r) +
                _symm_s(s) +
                _symm_s(t)
                )
            self.degree = 7
        elif index == 9:
            # Stroud number C2:7-3.
            r = numpy.sqrt(12.0 / 35.0)
            c = 3.0 * numpy.sqrt(186.0)
            s = numpy.sqrt((93.0 + c) / 155.0)
            t = numpy.sqrt((93.0 - c) / 155.0)
            w1 = 8.0 / 162.0
            w2 = 98.0 / 162.0
            w3 = 31.0 / 162.0
            self.weights = (
                [w1] +
                4 * [w2] +
                8 * [w3]
                )
            self.points = (
                [[0.0, 0.0]] +
                _symm_r_0(r) +
                _symm_s_t(s, t)
                )
            self.degree = 7
        else:
            assert index == 10
            scheme1d = line_segment.GaussLegendre(8)
            self.weights = numpy.outer(
                scheme1d.weights, scheme1d.weights
                ).flatten()
            self.points = numpy.dstack(numpy.meshgrid(
                scheme1d.points, scheme1d.points
                )).reshape(-1, 2)
            assert len(self.points) == 64
            self.degree = 15

        self.weights = numpy.array(self.weights)
        self.points = numpy.array(self.points)
        return
