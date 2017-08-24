# -*- coding: utf-8 -*-
#
import numpy

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
        elif index == 'Tn 2-1a':
            self.set_data(HammerStroud(n, '1a'))
        elif index == 'Tn 2-1b':
            self.set_data(HammerStroud(n, '1b'))
        elif index == 'Tn 2-2':
            self.set_data(Lauffer(n, 2))
        elif index == 'Tn 3-1':
            self.set_data(HammerStroud(n, '2'))
        elif index == 'Tn 3-2':
            self.set_data(Stroud1966(n, 'I'))
        elif index == 'Tn 3-3':
            self.set_data(Stroud1961(n))
        elif index == 'Tn 3-4':
            self.set_data(Stroud1966(n, 'II'))
        elif index == 'Tn 3-5':
            self.set_data(Stroud1966(n, 'III'))
        elif index == 'Tn 3-6a':
            self.set_data(Stroud1964(n, variant='a'))
        elif index == 'Tn 3-6b':
            self.set_data(Stroud1964(n, variant='b'))
        elif index == 'Tn 3-7':
            self.set_data(Stroud1966(n, 'IV'))
        elif index == 'Tn 3-8':
            self.set_data(Stroud1966(n, 'V'))
        elif index == 'Tn 3-9':
            self.set_data(Lauffer(n, 3))
        elif index == 'Tn 3-10':
            self.set_data(Stroud1966(n, 'VI'))
        elif index == 'Tn 3-11':
            self.set_data(Stroud1966(n, 'VII'))
        elif index == 'Tn 4-1':
            self.set_data(Lauffer(n, 4))
        elif index == 'Tn 5-1':
            self.set_data(Stroud1969(n))
        else:
            assert index == 'Tn 5-2'
            self.set_data(Lauffer(n, 5))
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.bary = scheme.bary
        self.points = scheme.points
        return
