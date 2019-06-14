# -*- coding: utf-8 -*-
#
"""
J. Albrecht, L. Collatz,
Zur numerischen Auswertung mehrdimensionaler Integrale,
ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
<https://doi.org/10.1002/zamm.19580380102>
"""
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, pm
from .helpers import DiskScheme


def AlbrechtCollatz(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi

    data = [
        # ERR Wrongly stated in Stroud as sqrt(1/2) instead of 1/2
        (frac(1, 4), pm(2, frac(1, 2)))
    ]
    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Albrecht-Collatz", 3, weights, points)
