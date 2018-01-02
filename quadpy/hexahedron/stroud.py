# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .albrecht_collatz import AlbrechtCollatz
from .hammer_stroud import HammerStroud
from .hammer_wymore import HammerWymore
from .mustard_lyness_blatt import MustardLynessBlatt
from .sadowsky import Sadowsky
from .sarma_stroud import SarmaStroud
from .stroud1967 import Stroud1967
from .tyler import Tyler

from .helpers import pm_rrr

from ..ncube import Ewing
from ..helpers import untangle


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        reference_volume = 8
        if index == 'C3 3-1':
            self.set_data(Tyler(1))
        elif index == 'C3 3-2':
            # product Gauss
            self.degree = 3
            data = [
                (frac(1, 8), pm_rrr(sqrt(frac(1, 3))))
                ]
            self.points, self.weights = untangle(data)
            self.weights *= reference_volume
        elif index == 'C3 3-3':
            self.set_data(Ewing(3, symbolic))
        elif index == 'C3 3-4':
            self.set_data(MustardLynessBlatt(1, symbolic))
        elif index == 'C3 3-5':
            self.set_data(MustardLynessBlatt(2, symbolic))
        elif index == 'C3 3-6':
            self.set_data(AlbrechtCollatz(symbolic))
        elif index == 'C3 3-7':
            self.set_data(MustardLynessBlatt(3, symbolic))
        elif index == 'C3 5-1':
            self.set_data(Stroud1967(symbolic))
        elif index == 'C3 5-2':
            self.set_data(HammerStroud('2-3', symbolic))
        elif index == 'C3 5-3':
            self.set_data(Tyler(2, symbolic))
        elif index == 'C3 5-4':
            self.set_data(MustardLynessBlatt(4, symbolic))
        elif index == 'C3 5-5':
            self.set_data(MustardLynessBlatt(5, symbolic))
        elif index == 'C3 5-6':
            self.set_data(MustardLynessBlatt(6, symbolic))
        elif index == 'C3 5-7':
            self.set_data(MustardLynessBlatt(7, symbolic))
        elif index == 'C3 5-8':
            self.set_data(Sadowsky(symbolic))
        elif index == 'C3 7-1a':
            self.set_data(HammerStroud('5-3a', symbolic))
        elif index == 'C3 7-1b':
            self.set_data(HammerStroud('5-3b', symbolic))
        elif index == 'C3 7-2':
            self.set_data(HammerWymore(symbolic=symbolic))
        else:
            assert index == 'C3 7-3'
            self.set_data(SarmaStroud(symbolic))
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
