# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

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
            (fr(250, 225), _z()),
            (-fr(8, 225), _symm_r_0(1)),
            (fr(7, 900), _symm_s(1)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
