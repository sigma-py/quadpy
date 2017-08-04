# -*- coding: utf-8 -*-
#
import warnings

import numpy

from .albrecht_collatz import AlbrechtCollatz
from .burnside import Burnside
from .irwin import Irwin
from .maxwell import Maxwell
from .meister import Meister
from .miller import Miller
from .phillips import Phillips
from .rabinowitz_richter import RabinowitzRichter
from .tyler import Tyler
from .helpers import _symm_r_0, _symm_s, _symm_s_t, _z

from .. import ncube
from ..helpers import untangle


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
            r = numpy.sqrt(3.0 / 5.0)
            data = [
                (16.0/81.0, _z()),
                (10.0/81.0, _symm_r_0(r)),
                (25.0/324.0, _symm_s(r)),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= reference_volume
        elif index == 'C2 5-5':
            self.set_data(Tyler(1))
        elif index == 'C2 5-6':
            self.set_data(AlbrechtCollatz(4))
        elif index == 'C2 5-7':
            self.set_data(Irwin(2))
        elif index == 'C2 7-1':
            self.set_data(Tyler(2))
        elif index == 'C2 7-2':
            self.set_data(Phillips())
        elif index == 'C2 7-3':
            self.set_data(Maxwell())
        elif index == 'C2 7-4':
            # product Gauss
            # TODO fix
            warnings.warn('Formula {} only has degree 1!'.format(index))
            self.degree = 1

            r = numpy.sqrt((15.0 - 2.0*numpy.sqrt(30)) / 35.0)
            s = numpy.sqrt((15.0 + 2.0*numpy.sqrt(30)) / 35.0)

            B1 = (59.0 + 6.0*numpy.sqrt(30.0)) / 864.0
            B2 = (59.0 - 6.0*numpy.sqrt(30.0)) / 864.0
            B3 = 49.0 / 864.0

            r = numpy.sqrt(3.0 / 5.0)
            data = [
                (B1, _symm_s(r)),
                (B2, _symm_s(s)),
                (B3, _symm_s_t(r, s)),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= reference_volume
        elif index == 'C2 7-5':
            self.set_data(Tyler(3))
        elif index == 'C2 7-6':
            self.set_data(Meister())
        elif index == 'C2 9-1':
            self.set_data(RabinowitzRichter(1))
        elif index == 'C2 11-1':
            self.set_data(RabinowitzRichter(2))
        elif index == 'C2 11-2':
            self.set_data(RabinowitzRichter(3))
        elif index == 'C2 13-1':
            self.set_data(RabinowitzRichter(4))
        elif index == 'C2 15-1':
            self.set_data(RabinowitzRichter(5))
        else:
            assert index == 'C2 15-2'
            self.set_data(RabinowitzRichter(6))

        # assert index == 'C2 5-7', 'Illegal index \'{}\'.'.format(index)

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
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
