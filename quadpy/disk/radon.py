# -*- coding: utf-8 -*-
#
"""
J. Radon,
Zur mechanischen Kubatur,
Monatshefte für Mathematik (1948),
Volume: 52, page 286-300, ISSN: 0026-9255; 1436-5081/e.
<https://eudml.org/doc/176796>.
"""
from __future__ import division

import numpy
import sympy

from ..helpers import z, untangle
from .helpers import DiskScheme


def Radon(alpha, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    r = sqrt(frac(alpha + 4, alpha + 6))
    s = sqrt(frac(alpha + 4, 4 * (alpha + 6)))
    t = sqrt(frac(3 * (alpha + 4), 4 * (alpha + 6)))

    A = frac(4, (alpha + 4) ** 2)
    B = frac((alpha + 2) * (alpha + 6), 6 * (alpha + 4) ** 2)

    data = [
        (A, z(2)),
        (B, numpy.array([[+r, 0], [-r, 0]])),
        # Stroud is missing +- in front of t.
        (B, numpy.array([[+s, +t], [-s, +t], [+s, -t], [-s, -t]])),
    ]

    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Radon", 5, weights, points)
