# -*- coding: utf-8 -*-
#
import numpy
import sympy

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
    def __init__(self, index, symbolic=True):
        self.name = 'Stroud_U3({})'.format(index)

        if index == 'U3 3-1':
            self.set_data(McLaren(1, symbolic=symbolic))
        elif index == 'U3 5-1':
            self.set_data(AlbrechtCollatz(1, symbolic=symbolic))
        elif index == 'U3 5-2':
            self.set_data(AlbrechtCollatz(2, symbolic=symbolic))
        elif index == 'U3 5-3':
            self.set_data(AlbrechtCollatz(3, symbolic=symbolic))
        elif index == 'U3 5-4':
            self.set_data(AlbrechtCollatz(4, symbolic=symbolic))
        elif index == 'U3 5-5':
            self.set_data(McLaren(2, symbolic=symbolic))
        elif index == 'U3 7-1':
            self.set_data(McLaren(3, symbolic=symbolic))
        elif index == 'U3 7-2':
            self.set_data(AlbrechtCollatz(5, symbolic=symbolic))
        elif index == 'U3 8-1':
            self.set_data(McLaren(4, symbolic=symbolic))
        elif index == 'U3 9-1':
            self.set_data(McLaren(5, symbolic=symbolic))
        elif index == 'U3 9-2':
            self.set_data(McLaren(6, symbolic=symbolic))
        elif index == 'U3 9-3':
            self.set_data(McLaren(7, symbolic=symbolic))
        elif index == 'U3 11-1':
            self.set_data(McLaren(8, symbolic=symbolic))
        elif index == 'U3 11-2':
            scheme = Stroud1969(3, symbolic=symbolic)
            self.degree = scheme.degree
            self.weights = scheme.weights
            pi = sympy.pi if symbolic else numpy.pi
            self.weights /= 4 * pi
            self.points = scheme.points
            self.azimuthal_polar = cartesian_to_spherical(self.points)
        elif index == 'U3 11-3':
            self.set_data(McLaren(9, symbolic=symbolic))
        else:
            assert index == 'U3 14-1', 'Illegal index {}.'.format(index)
            self.set_data(McLaren(10, symbolic=symbolic))
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        self.azimuthal_polar = scheme.azimuthal_polar
        return
