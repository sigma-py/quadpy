# -*- coding: utf-8 -*-
#
from .stroud1967 import Stroud1967
from .stroud_secrest import StroudSecrest


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, index):
        self.name = 'Stroud_Enr2({})'.format(index)
        self.dim = n
        if index == '3-1':
            self.set_data(StroudSecrest(n, 'I'))
        elif index == '3-2':
            self.set_data(StroudSecrest(n, 'III'))
        elif index == '5-1a':
            self.set_data(Stroud1967(n, 'a'))
        elif index == '5-1b':
            self.set_data(Stroud1967(n, 'b'))
        elif index == '5-2':
            self.set_data(StroudSecrest(n, 'IV'))
        elif index == '5-3':
            # TODO
            pass
        elif index == '5-4':
            # spherical produce Lobatto
            # TODO
            pass
        elif index == '5-5':
            # TODO
            pass
        elif index == '5-6':
            # TODO
            pass
        elif index == '7-1':
            # TODO
            pass
        elif index == '7-2':
            # TODO
            pass
        elif index == '7-3':
            # TODO
            pass
        elif index == '9-1':
            # TODO
            pass
        else:
            assert index == '11-1'
            # TODO
            pass
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
