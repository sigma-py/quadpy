# -*- coding: utf-8 -*-
#
from . import albrecht_collatz
from . import hammer_stroud
from . import mysovskih


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
            self.set_data(hammer_stroud.HammerStroud())
        elif index == 'S2 3-2':
            self.set_data(albrecht_collatz.AlbrechtCollatz())
        elif index == 'S2 4-1':
            self.set_data(mysovskih.Mysovskih(0.0))
        else:
            assert False

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
