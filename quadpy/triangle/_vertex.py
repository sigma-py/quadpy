# -*- coding: utf-8 -*-
#

import sympy

from ._helpers import s2, TriangleScheme


def vertex(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = s2([frac(1, 3), 0])
    return TriangleScheme("Vertex scheme", weights, points, 1)
