# -*- coding: utf-8 -*-
#
from . import stroud1957


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
        if index == 'Sn 2-1':
            scheme = stroud1957.Stroud1957(n)
            self.degree = scheme.degree
            self.points = scheme.points
            self.weights = scheme.weights
        else:
            assert False

        return
