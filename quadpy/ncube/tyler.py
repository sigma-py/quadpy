# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _z, _fsd


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://dx.doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self, n):
        reference_volume = 2.0**n
        self.degree = 3
        self.weights = numpy.concatenate([
            numpy.full(1, (3.0 - n)/3.0 * reference_volume),
            numpy.full(2*n, reference_volume/6.0),
            ])
        self.points = numpy.concatenate([
            _z(n),
            _fsd(n, 1.0, 1)
            ])
        return
