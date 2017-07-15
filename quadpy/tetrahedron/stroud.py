# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s4, _s31


class Stroud(object):
    '''
    A. H. Stroud,
    Approximate calculation of multiple integrals,
    Englewood Cliffs, NJ : Prentice-Hall, c 1971. - XIII,
    ISBN 0-13-043893-6.
    '''
    def __init__(self, index):
        self.name = 'Stroud({})'.format(index)
        if index == 0:
            self.degree = 2
            self.weights = numpy.full(4, 1.0/4.0)
            bary = _s31((5.0 - numpy.sqrt(5.0)) / 20.0)
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
