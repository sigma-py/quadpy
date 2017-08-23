# -*- coding: utf-8 -*-
#
import numpy

from .lauffer import Lauffer

from ..helpers import untangle


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
        if index == 'Tn 1-1':
            # midpoint formula
            self.degree = 1
            data = [
                (1.0, numpy.array([numpy.full(n+1, 1.0/(n+1))]))
                ]
            self.bary, self.weights = untangle(data)
            self.points = self.bary[:, 1:]
        elif index == 'Tn 1-2':
            self.set_data(Lauffer(n, 1))
        else:
            assert False
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
