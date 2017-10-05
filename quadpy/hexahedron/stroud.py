# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

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
    def __init__(self, index):
        reference_volume = 8
        if index == 'C3 3-1':
            self.set_data(Tyler(1))
        elif index == 'C3 3-2':
            # product Gauss
            self.degree = 3
            data = [
                (fr(1, 8), pm_rrr(sqrt(fr(1, 3))))
                ]
            self.points, self.weights = untangle(data)
            self.weights *= reference_volume
        elif index == 'C3 3-3':
            self.set_data(Ewing(3))
        elif index == 'C3 3-4':
            self.set_data(MustardLynessBlatt(1))
        elif index == 'C3 3-5':
            self.set_data(MustardLynessBlatt(2))
        elif index == 'C3 3-6':
            self.set_data(AlbrechtCollatz())
        elif index == 'C3 3-7':
            self.set_data(MustardLynessBlatt(3))
        elif index == 'C3 5-1':
            self.set_data(Stroud1967())
        elif index == 'C3 5-2':
            self.set_data(HammerStroud('2-3'))
        elif index == 'C3 5-3':
            self.set_data(Tyler(2))
        elif index == 'C3 5-4':
            self.set_data(MustardLynessBlatt(4))
        elif index == 'C3 5-5':
            self.set_data(MustardLynessBlatt(5))
        elif index == 'C3 5-6':
            self.set_data(MustardLynessBlatt(6))
        elif index == 'C3 5-7':
            self.set_data(MustardLynessBlatt(7))
        elif index == 'C3 5-8':
            self.set_data(Sadowsky())
        elif index == 'C3 7-1a':
            self.set_data(HammerStroud('5-3a'))
        elif index == 'C3 7-1b':
            self.set_data(HammerStroud('5-3b'))
        elif index == 'C3 7-2':
            self.set_data(HammerWymore())
        else:
            assert index == 'C3 7-3'
            self.set_data(SarmaStroud())
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
