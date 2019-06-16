# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import s2, TriangleScheme


def Vertex(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = s2([frac(1, 3), 0])
    points = bary[:, 1:]
    return TriangleScheme("Vertex scheme", 1, weights, points, bary)
