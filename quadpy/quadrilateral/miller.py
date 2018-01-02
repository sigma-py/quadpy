# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import _symm_r_0, _symm_s, _z
from ..helpers import untangle


class Miller(object):
    '''
    J.C.P. Miller,
    Numerical Quadrature Over a Rectangular Domain in Two or More Dimensions.
    Part 3: Quadrature of a Harmonic Integrand,
    Mathematics of Computation,
    Vol. 14, No. 71 (Jul., 1960), pp. 240-248,
    <https://doi.org/10.2307/2003163>.

    This scheme is exact for harmonic integrands of degree <= 11.
    '''
    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y

        self.name = 'Miller'
        self.degree = 1
        data = [
            (frac(250, 225), _z()),
            (-frac(8, 225), _symm_r_0(1)),
            (frac(7, 900), _symm_s(1)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return
