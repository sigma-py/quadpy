# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _symm_r_0, _symm_s, _z


class Miller(object):
    '''
    J.C.P. Miller,
    Numerical Quadrature Over a Rectangular Domain in Two or More Dimensions.
    Part 3: Quadrature of a Harmonic Integrand,
    Mathematics of Computation,
    Vol. 14, No. 71 (Jul., 1960), pp. 240-248,
    <https://dx.doi.org/10.2307/2003163>.
    '''
    def __init__(self):
        self.name = 'Miller'
        self.degree = 1
        self.weights = numpy.concatenate([
            numpy.full(1, 250.0/225.0),
            numpy.full(4, -8.0/225.0),
            numpy.full(4, 7.0/900.0),
            ])
        self.points = numpy.concatenate([
            _z(),
            _symm_r_0(1.0),
            _symm_s(1.0)
            ])
        return
