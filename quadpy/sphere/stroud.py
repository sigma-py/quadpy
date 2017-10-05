# -*- coding: utf-8 -*-
#
import numpy

from .helpers import cartesian_to_spherical

from .albrecht_collatz import AlbrechtCollatz
from .mclaren import McLaren
from ..nsphere.stroud1969 import Stroud1969


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        if index == 'U3 3-1':
            self.set_data(McLaren(1))
        elif index == 'U3 5-1':
            self.set_data(AlbrechtCollatz(1))
        elif index == 'U3 5-2':
            self.set_data(AlbrechtCollatz(2))
        elif index == 'U3 5-3':
            self.set_data(AlbrechtCollatz(3))
        elif index == 'U3 5-4':
            self.set_data(AlbrechtCollatz(4))
        elif index == 'U3 5-5':
            self.set_data(McLaren(2))
        elif index == 'U3 7-1':
            self.set_data(McLaren(3))
        elif index == 'U3 7-2':
            self.set_data(AlbrechtCollatz(5))
        elif index == 'U3 8-1':
            self.set_data(McLaren(4))
        elif index == 'U3 9-1':
            self.set_data(McLaren(5))
        elif index == 'U3 9-2':
            self.set_data(McLaren(6))
        elif index == 'U3 9-3':
            self.set_data(McLaren(7))
        elif index == 'U3 11-1':
            self.set_data(McLaren(8))
        elif index == 'U3 11-2':
            scheme = Stroud1969(3)
            self.degree = scheme.degree
            self.weights = scheme.weights
            self.weights /= 4 * numpy.pi
            flt = numpy.vectorize(float)
            self.points = flt(scheme.points)
            self.phi_theta = cartesian_to_spherical(self.points)
        elif index == 'U3 11-3':
            self.set_data(McLaren(9))
        else:
            assert index == 'U3 14-1', 'Illegal index {}.'.format(index)
            self.set_data(McLaren(10))
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        self.phi_theta = scheme.phi_theta
        return
