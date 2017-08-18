# -*- coding: utf-8 -*-
#
from .helpers import _symm_r_0, _symm_s, _z

from ..helpers import untangle


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
        # This scheme is exact for harmonic integrands of degree <= 11.
        self.name = 'Miller'
        self.degree = 1
        data = [
            (250.0/225.0, _z()),
            (-8.0/225.0, _symm_r_0(1.0)),
            (7.0/900.0, _symm_s(1.0)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4.0
        return
