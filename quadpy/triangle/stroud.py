# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from __future__ import division

import numpy

from .albrecht_collatz import AlbrechtCollatz
from .hammer_marlowe_stroud import HammerMarloweStroud

from .helpers import TriangleScheme
from ..line_segment.gauss_legendre import GaussLegendre


def stroud_7_1():
    # conical product Gauss
    gl4 = GaussLegendre(4)
    r = (gl4.points + 1) / 2
    A = gl4.weights / 2

    # Generate Gauss formula for int_0^1 (1-s) * f(s) ds.
    # ```
    # k = numpy.arange(8)
    # moments = 1 / (k**2 + 3*k + 2)
    # alpha, beta = orthopy.line.chebyshev(moments)
    # s, B = orthopy.line.schemes.custom(alpha, beta, mode='numpy')
    # ```
    s = numpy.array(
        [
            5.710419611452533e-02,
            2.768430136381415e-01,
            5.835904323689318e-01,
            8.602401356562251e-01,
        ]
    )
    B = numpy.array(
        [
            1.355069134315012e-01,
            2.034645680102685e-01,
            1.298475476082247e-01,
            3.118097095000554e-02,
        ]
    )

    weights = numpy.array([2 * A[i] * B[j] for i in range(4) for j in range(4)])
    bary = numpy.array(
        [
            [s[j], r[i] * (1 - s[j]), (1 - r[i]) * (1 - s[j])]
            for i in range(4)
            for j in range(4)
        ]
    )
    return TriangleScheme("Stroud 7-1", 7, weights, bary)


Stroud = {
    "T2 3-1": AlbrechtCollatz,
    "T2 5-1": HammerMarloweStroud[5],
    "T2 7-1": stroud_7_1,
}
