# -*- coding: utf-8 -*-
#
from .stroud_secrest import StroudSecrest


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, index):
        self.name = 'Stroud({})'.format(index)
        self.dim = n
        if index == 'Enr2 3-1':
            self.set_data(StroudSecrest('I', n))
        else:
            assert False
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
