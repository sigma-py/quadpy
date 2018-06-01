# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

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
    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        pi = sympy.pi if symbolic else numpy.pi
        cos = numpy.vectorize(sympy.cos) if symbolic else numpy.cos
        sin = numpy.vectorize(sympy.sin) if symbolic else numpy.sin
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        pm_ = numpy.array([+1, -1])

        self.name = 'Stroud({})'.format(index)
        if index == 'S2 3-1':
            self.set_data(
                hammer_stroud.HammerStroud('11-2', symbolic=symbolic)
                )
        elif index == 'S2 3-2':
            self.set_data(albrecht_collatz.AlbrechtCollatz(symbolic=symbolic))
        elif index == 'S2 4-1':
            self.set_data(mysovskih.Mysovskih(1, symbolic=symbolic))
        elif index == 'S2 5-1':
            self.set_data(radon.Radon(0, symbolic=symbolic))
        elif index == 'S2 5-2':
            self.degree = 5
            r = sqrt(frac(1, 2))
            data = [
                (frac(1, 6), z(2)),
                (frac(1, 6), fsd(2, (r, 1))),
                (frac(1, 24), pm(2, r))
                ]
            self.points, self.weights = untangle(data)
            self.weights *= pi
        elif index == 'S2 7-1':
            self.set_data(peirce1956.Peirce1956(1, symbolic=symbolic))
        elif index == 'S2 7-2':
            # spherical product Gauss
            self.degree = 7

            r1, r2 = sqrt((3 - pm_ * sqrt(3)) / 6)

            a = (2*numpy.arange(8)+1) * pi / 8
            x = numpy.array([cos(a), sin(a)]).T

            data = [
                (frac(1, 16), r1*x),
                (frac(1, 16), r2*x),
                ]
            self.points, self.weights = untangle(data)
            self.weights *= pi
        elif index == 'S2 9-1':
            self.set_data(albrecht.Albrecht(4, symbolic=symbolic))
        elif index == 'S2 9-2':
            self.set_data(rabinowitz_richter.RabinowitzRichter(1))
        elif index == 'S2 9-3':
            # spherical product Gauss
            self.degree = 9

            r1, r2 = sqrt((6 - pm_ * sqrt(6)) / 10)

            a = (numpy.arange(10)+1) * pi / 5
            x = numpy.array([cos(a), sin(a)]).T

            B0 = frac(1, 9)
            B1, B2 = (16 + pm_ * sqrt(6)) / 360

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
            self.set_data(peirce1956.Peirce1956(2, symbolic=symbolic))
        elif index == 'S2 11-1':
            self.set_data(mysovskih.Mysovskih(2, symbolic=symbolic))
        elif index == 'S2 11-2':
            self.set_data(albrecht.Albrecht(5, symbolic=symbolic))
        elif index == 'S2 11-3':
            self.set_data(rabinowitz_richter.RabinowitzRichter(4))
        elif index == 'S2 11-4':
            self.set_data(peirce1956.Peirce1956(3, symbolic=symbolic))
        elif index == 'S2 13-1':
            self.set_data(rabinowitz_richter.RabinowitzRichter(5))
        elif index == 'S2 13-2':
            self.set_data(albrecht.Albrecht(6, symbolic=symbolic))
        elif index == 'S2 15-1':
            self.set_data(mysovskih.Mysovskih(3, symbolic=symbolic))
        elif index == 'S2 15-2':
            self.set_data(albrecht.Albrecht(7, symbolic=symbolic))
        else:
            assert index == 'S2 17-1'
            self.set_data(albrecht.Albrecht(8, symbolic=symbolic))

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
