# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _z, _pm


class Ewing(object):
    '''
    G.M. Ewing,
    On Approximate Cubature,
    The American Mathematical Monthly,
    Vol. 48, No. 2 (Feb., 1941), pp. 134-136,
    <https://dx.doi.org/dx.doi.org/10.2307/2303604>.
    '''
    def __init__(self, n):
        reference_volume = 2.0**n
        self.degree = 3
        self.weights = numpy.concatenate([
            numpy.full(1, 2.0/3.0 * reference_volume),
            numpy.full(2**n, 1.0/3.0 / 2**n * reference_volume),
            ])
        self.points = numpy.concatenate([
            _z(n),
            _pm(n, 1.0),
            ])
        return
