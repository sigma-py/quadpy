# -*- coding: utf-8 -*-
#
"""
J. Albrecht, L. Collatz,
Zur numerischen Auswertung mehrdimensionaler Integrale,
ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
<https://doi.org/10.1002/zamm.19580380102>
"""
from __future__ import division

import sympy

from .helpers import untangle2, TriangleScheme


def AlbrechtCollatz(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = {"s2": [[frac(2, 30), frac(1, 2)], [frac(9, 15), frac(1, 6)]]}
    bary, weights = untangle2(data)

    points = bary[:, 1:]
    weights /= 2
    return TriangleScheme("Albrecht-Collatz", 3, weights, points, bary)
