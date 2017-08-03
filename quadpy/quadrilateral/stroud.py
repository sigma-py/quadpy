# -*- coding: utf-8 -*-
#
import numpy

from .albrecht_collatz import AlbrechtCollatz
from .burnside import Burnside
from .irwin import Irwin
from .miller import Miller
from .tyler import Tyler
from .helpers import _symm_r_0, _symm_s, _symm_s_t, _z

from .. import line_segment
from .. import ncube


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.

    <https://people.sc.fsu.edu/~jburkardt/m_src/stroud/square_unit_set.m>
    <http://nines.cs.kuleuven.be/ecf/mtables.html>
    '''
    def __init__(self, index):
        reference_volume = 4.0
        if index == 'C2 1-1':
            # Product trapezoidal
            self.degree = 1
            self.weights = reference_volume * numpy.full(4, 0.25)
            self.points = _symm_s(1.0)
        elif index == 'C2 1-2':
            self.set_data(Miller())
        elif index == 'C2 3-1':
            # Product Gauss
            self.degree = 3
            self.weights = reference_volume * numpy.full(4, 0.25)
            # misprint in Stroud: sqrt(1/3) vs 1/3
            self.points = _symm_s(numpy.sqrt(1.0/3.0))
        elif index == 'C2 3-2':
            self.set_data(ncube.Ewing(2))
        elif index == 'C2 3-3':
            # product Simpson
            self.set_data(ncube.Stroud(2, 'Cn 3-6'))
        elif index == 'C2 3-4':
            self.set_data(AlbrechtCollatz(1))
        elif index == 'C2 3-5':
            self.set_data(Irwin(1))
        elif index == 'C2 5-1':
            self.set_data(AlbrechtCollatz(2))
        elif index == 'C2 5-2':
            self.set_data(AlbrechtCollatz(3))
        elif index == 'C2 5-3':
            self.set_data(Burnside())
        elif index == 'C2 5-4':
            # product Gauss
            self.degree = 5
            self.weights = reference_volume * numpy.concatenate([
                numpy.full(1, 16.0/81.0),
                numpy.full(4, 10.0/81.0),
                numpy.full(4, 25.0/324.0)
                ])
            r = numpy.sqrt(3.0 / 5.0)
            self.points = numpy.concatenate([
                _z(),
                _symm_r_0(r),
                _symm_s(r)
                ])
        elif index == 'C2 5-5':
            self.set_data(Tyler())
        elif index == 'C2 5-6':
            self.set_data(AlbrechtCollatz(4))
        else:
            assert index == 'C2 5-7', 'Illegal index \'{}\'.'.format(index)
            self.set_data(Irwin(2))

        # elif index == 2:
        #     self.weights = 4 * [1.0]
        #     self.points = _symm_r_0(numpy.sqrt(2.0/3.0))
        #     self.degree = 3
        # elif index == 3:
        #     self.weights = (
        #         4 * [2.0/3.0] +
        #         [4.0/3.0]
        #         )
        #     self.points = (
        #         _symm_r_0(1.0) +
        #         [[0.0, 0.0]]
        #         )
        #     self.degree = 3
        # elif index == 4:
        #     self.weights = (
        #         4 * [1.0/3.0] +
        #         [8.0/3.0]
        #         )
        #     self.points = (
        #         _symm_s(1.0) +
        #         [[0.0, 0.0]]
        #         )
        #     self.degree = 3
        # elif index == 5:
        #     self.weights = (
        #         4 * [40.0 / 49.0] +
        #         4 * [9.0 / 49.0]
        #         )
        #     self.points = (
        #         _symm_r_0(numpy.sqrt(7.0/15.0)) +
        #         _symm_s(numpy.sqrt(7.0/9.0))
        #         )
        #     self.degree = 5
        # elif index == 6:
        #     self.weights = (
        #         4 * [1.0 / 9.0] +
        #         1 * [-8.0 / 9.0] +
        #         4 * [10.0 / 9.0]
        #         )
        #     self.points = (
        #         _symm_s(1.0) +
        #         [[0.0, 0.0]] +
        #         _symm_r_0(numpy.sqrt(2.0 / 5.0))
        #         )
        #     self.degree = 5
        # elif index == 7:
        #     self.weights = (
        #         [8.0 / 7.0] +
        #         8 * [5.0 / 14.0]
        #         )
        #     self.points = (
        #         [[0.0, 0.0]] +
        #         _symm_s_t(
        #             0.846233119448533574334773553209421,
        #             0.466071712166418914664132375714293
        #             )
        #         )
        #     self.degree = 5
        # elif index == 8:
        #     # Stroud number C2:7-1.
        #     r = numpy.sqrt(6.0 / 7.0)
        #     c = 3.0 * numpy.sqrt(583.0)
        #     s = numpy.sqrt((114.0 - c) / 287.0)
        #     t = numpy.sqrt((114.0 + c) / 287.0)
        #     w1 = 4.0 * 49.0 / 810.0
        #     w2 = 4.0 * (178981.0 + 923.0 * c) / 1888920.0
        #     w3 = 4.0 * (178981.0 - 923.0 * c) / 1888920.0
        #     #
        #     self.weights = (
        #         4 * [w1] +
        #         4 * [w2] +
        #         4 * [w3]
        #         )
        #     self.points = (
        #         _symm_r_0(r) +
        #         _symm_s(s) +
        #         _symm_s(t)
        #         )
        #     self.degree = 7
        # elif index == 9:
        #     # Stroud number C2:7-3.
        #     r = numpy.sqrt(12.0 / 35.0)
        #     c = 3.0 * numpy.sqrt(186.0)
        #     s = numpy.sqrt((93.0 + c) / 155.0)
        #     t = numpy.sqrt((93.0 - c) / 155.0)
        #     w1 = 8.0 / 162.0
        #     w2 = 98.0 / 162.0
        #     w3 = 31.0 / 162.0
        #     self.weights = (
        #         [w1] +
        #         4 * [w2] +
        #         8 * [w3]
        #         )
        #     self.points = (
        #         [[0.0, 0.0]] +
        #         _symm_r_0(r) +
        #         _symm_s_t(s, t)
        #         )
        #     self.degree = 7
        # else:
        #     assert index == 10
        #     scheme1d = line_segment.GaussLegendre(8)
        #     self.weights = numpy.outer(
        #         scheme1d.weights, scheme1d.weights
        #         ).flatten()
        #     self.points = numpy.dstack(numpy.meshgrid(
        #         scheme1d.points, scheme1d.points
        #         )).reshape(-1, 2)
        #     assert len(self.points) == 64
        #     self.degree = 15

        self.weights = numpy.array(self.weights)
        self.points = numpy.array(self.points)
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
