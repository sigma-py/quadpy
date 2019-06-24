# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from ._helpers import s2, TriangleScheme


def vertex(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = s2([frac(1, 3), 0])
    return TriangleScheme("Vertex scheme", weights, bary, 1)
