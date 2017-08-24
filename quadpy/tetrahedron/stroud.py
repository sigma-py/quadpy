# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .helpers import _s4, _s31

from ..nsimplex.stroud import Stroud as nsimplex

class Stroud(object):
    '''
    A. H. Stroud,
    Approximate calculation of multiple integrals,
    Englewood Cliffs, NJ : Prentice-Hall, c 1971. - XIII,
    ISBN 0-13-043893-6.
    '''
    def __init__(self, index):
        self.name = 'Stroud({})'.format(index)
        if index == 'T3 5-1':
            self.degree = 5
        else:
            assert index == 1
            self.degree = 3
            self.weights = 6 * numpy.concatenate([
                numpy.full(1, -2.0/15.0),
                numpy.full(4, 3.0/40.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/6.0)
                ])

        self.points = bary[:, 1:]
        return
