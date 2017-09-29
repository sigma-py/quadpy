# -*- coding: utf-8 -*-
#
import numpy
from sympy import pi, Rational as fr, sqrt, cos, sin

from . import albrecht
from . import albrecht_collatz
from . import hammer_stroud
from . import mysovskih
from . import peirce1956
from . import rabinowitz_richter
from . import radon

from ..helpers import z, fsd, pm, untangle


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        self.name = 'Stroud({})'.format(index)
        if index == 'S2 3-1':
            self.set_data(hammer_stroud.HammerStroud('11-2'))
        elif index == 'S2 3-2':
            self.set_data(albrecht_collatz.AlbrechtCollatz())
        elif index == 'S2 4-1':
            self.set_data(mysovskih.Mysovskih(1))
        elif index == 'S2 5-1':
            self.set_data(radon.Radon(0.0))
        elif index == 'S2 5-2':
            self.degree = 5
            r = sqrt(fr(1, 2))
            data = [
                (fr(1, 6), z(2)),
                (fr(1, 6), fsd(2, (r, 1))),
                (fr(1, 24), pm(2, r))
                ]
            self.points, self.weights = untangle(data)
            self.weights *= pi
        elif index == 'S2 7-1':
            self.set_data(peirce1956.Peirce1956(1))
        elif index == 'S2 7-2':
            # spherical product Gauss
            self.degree = 7

            r1, r2 = [sqrt((3 - t*sqrt(3)) / 6) for t in [+1, -1]]

            x = numpy.column_stack([
                [cos((2*k-1)*pi/8) for k in range(1, 9)],
                [sin((2*k-1)*pi/8) for k in range(1, 9)],
                ])

            data = [
                (fr(1, 16), r1*x),
                (fr(1, 16), r2*x),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= pi
        elif index == 'S2 9-1':
            self.set_data(albrecht.Albrecht(4))
        elif index == 'S2 9-2':
            self.set_data(rabinowitz_richter.RabinowitzRichter(1))
        elif index == 'S2 9-3':
            # spherical product Gauss
            self.degree = 9

            r1, r2 = [sqrt((6 - t*sqrt(6)) / 10) for t in [+1, -1]]

            x = numpy.column_stack([
                [cos(k * pi / 5) for k in range(1, 11)],
                [sin(k * pi / 5) for k in range(1, 11)],
                ])

            B0 = fr(1, 9)
            B1, B2 = [(16 + t*sqrt(6)) / 360 for t in [+1, -1]]

            data = [
                (B0, z(2)),
                (B1, r1*x),
                (B2, r2*x),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= pi
        elif index == 'S2 9-4':
            self.set_data(rabinowitz_richter.RabinowitzRichter(2))
        elif index == 'S2 9-5':
            self.set_data(peirce1956.Peirce1956(2))
        elif index == 'S2 11-1':
            self.set_data(mysovskih.Mysovskih(2))
        elif index == 'S2 11-2':
            self.set_data(albrecht.Albrecht(5))
        elif index == 'S2 11-3':
            self.set_data(rabinowitz_richter.RabinowitzRichter(4))
        elif index == 'S2 11-4':
            self.set_data(peirce1956.Peirce1956(3))
        elif index == 'S2 13-1':
            self.set_data(rabinowitz_richter.RabinowitzRichter(5))
        elif index == 'S2 13-2':
            self.set_data(albrecht.Albrecht(6))
        elif index == 'S2 15-1':
            self.set_data(mysovskih.Mysovskih(3))
        elif index == 'S2 15-2':
            self.set_data(albrecht.Albrecht(7))
        else:
            assert index == 'S2 17-1'
            self.set_data(albrecht.Albrecht(8))

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
