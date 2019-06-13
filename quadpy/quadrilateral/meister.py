# -*- coding: utf-8 -*-
#
"""
Bernd Meister,
On a Family of Cubature Formulae,
Comput J (1966) 8 (4): 368-371,
<https://doi.org/10.1093/comjnl/8.4.368>.
"""
from __future__ import division

import sympy

from .helpers import concat, zero, symm_s, symm_r0, symm_s_t, QuadrilateralScheme


def Meister(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    name = "Meister"
    degree = 7

    r = frac(2, 3)
    s = frac(1, 3)

    weights, points = concat(
        zero(frac(1024, 6720)),
        symm_s([frac(576, 6720), r], [-frac(9, 6720), s], [frac(47, 6720), 1]),
        symm_r0([frac(576, 6720), r]),
        symm_s_t([frac(117, 6720), 1, s]),
    )

    weights *= 4
    return QuadrilateralScheme(name, degree, weights, points)
