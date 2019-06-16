# -*- coding: utf-8 -*-
#
"""
Michael Sadowsky,
A Formula for Approximate Computation of a Triple Integral,
The American Mathematical Monthly,
Vol. 47, No. 8 (Oct., 1940), pp. 539-543,
<https://doi.org/10.2307/2303834>.
"""
from __future__ import division

import numpy
import sympy

from .helpers import fs_r00, fs_rr0, fs_rrs, HexahedronScheme
from ..helpers import untangle


def Sadowsky(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    data = [
        (frac(91, 450), fs_r00(1)),
        (frac(-20, 225), fs_rr0(1)),
        (frac(8, 225), fs_rrs(sqrt(frac(5, 8)), 1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Sadowsky", 5, weights, points)
