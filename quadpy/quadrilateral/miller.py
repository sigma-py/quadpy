# -*- coding: utf-8 -*-
#
"""
J.C.P. Miller,
Numerical Quadrature Over a Rectangular Domain in Two or More Dimensions.
Part 3: Quadrature of a Harmonic Integrand,
Mathematics of Computation,
Vol. 14, No. 71 (Jul., 1960), pp. 240-248,
<https://doi.org/10.2307/2003163>.

This scheme is exact for harmonic integrands of degree <= 11.
"""
from __future__ import division

import sympy

from .helpers import concat, zero, symm_r0, symm_s, QuadrilateralScheme


def Miller(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    name = "Miller"
    degree = 1

    weights, points = concat(
        zero(frac(250, 225)), symm_r0([-frac(8, 225), 1]), symm_s([frac(7, 900), 1])
    )
    weights *= 4
    return QuadrilateralScheme(name, degree, weights, points)
