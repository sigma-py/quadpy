# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .hammer_stroud import HammerStroud
from .lauffer import Lauffer
from .stroud1961 import Stroud1961
from .stroud1964 import Stroud1964
from .stroud1966 import Stroud1966
from .stroud1969 import Stroud1969

from ..helpers import untangle


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y

        self.name = 'Stroud({})'.format(index)
        self.dim = n
        if index == 'Tn 1-1':
            # midpoint formula
            self.degree = 1
            data = [
                (1, numpy.full((1, n+1), frac(1, n+1)))
                ]
            self.bary, self.weights = untangle(data)
            self.points = self.bary[:, 1:]
        elif index == 'Tn 1-2':
            self.set_data(Lauffer(n, 1, symbolic=symbolic))
        elif index == 'Tn 2-1a':
            self.set_data(HammerStroud(n, '1a', symbolic=symbolic))
        elif index == 'Tn 2-1b':
            self.set_data(HammerStroud(n, '1b', symbolic=symbolic))
        elif index == 'Tn 2-2':
            self.set_data(Lauffer(n, 2, symbolic=symbolic))
        elif index == 'Tn 3-1':
            self.set_data(HammerStroud(n, '2', symbolic=symbolic))
        elif index == 'Tn 3-2':
            self.set_data(Stroud1966(n, 'I', symbolic=symbolic))
        elif index == 'Tn 3-3':
            self.set_data(Stroud1961(n, symbolic=symbolic))
        elif index == 'Tn 3-4':
            self.set_data(Stroud1966(n, 'II', symbolic=symbolic))
        elif index == 'Tn 3-5':
            self.set_data(Stroud1966(n, 'III', symbolic=symbolic))
        elif index == 'Tn 3-6a':
            self.set_data(Stroud1964(n, variant='a', symbolic=symbolic))
        elif index == 'Tn 3-6b':
            self.set_data(Stroud1964(n, variant='b', symbolic=symbolic))
        elif index == 'Tn 3-7':
            self.set_data(Stroud1966(n, 'IV', symbolic=symbolic))
        elif index == 'Tn 3-8':
            self.set_data(Stroud1966(n, 'V', symbolic=symbolic))
        elif index == 'Tn 3-9':
            self.set_data(Lauffer(n, 3, symbolic=symbolic))
        elif index == 'Tn 3-10':
            self.set_data(Stroud1966(n, 'VI', symbolic=symbolic))
        elif index == 'Tn 3-11':
            self.set_data(Stroud1966(n, 'VII', symbolic=symbolic))
        elif index == 'Tn 4-1':
            self.set_data(Lauffer(n, 4, symbolic=symbolic))
        elif index == 'Tn 5-1':
            self.set_data(Stroud1969(n, symbolic=symbolic))
        else:
            assert index == 'Tn 5-2'
            self.set_data(Lauffer(n, 5, symbolic=symbolic))
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.bary = scheme.bary
        self.points = scheme.points
        return
