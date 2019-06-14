# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from __future__ import division

import numpy
import sympy

from .rabinowitz_richter import RabinowitzRichter
from .stroud_secrest import StroudSecrest
from .helpers import E2rScheme
from ..helpers import untangle


def stroud_4_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    cos = numpy.vectorize(sympy.cos) if symbolic else numpy.cos
    sin = numpy.vectorize(sympy.sin) if symbolic else numpy.sin
    pi = sympy.pi if symbolic else numpy.pi

    pts = (
        2
        * sqrt(5)
        * numpy.array(
            [
                [cos(2 * i * pi / 5) for i in range(5)],
                [sin(2 * i * pi / 5) for i in range(5)],
            ]
        ).T
    )
    data = [(frac(7, 10), numpy.array([[0, 0]])), (frac(3, 50), pts)]

    points, weights = untangle(data)
    weights *= 2 * pi
    return E2rScheme("Stroud 4-1", 4, weights, points)


# # The boolean tells whether the factor 2*pi is already in the weights
# _gen = {
#     "4-1": (_gen4_1, False),
#     "5-1": (stroud_secrest.v, False),
#     "7-1": (stroud_secrest.vi, False),
#     "9-1": (rabinowitz_richter.gen1, True),
#     "11-1": (rabinowitz_richter.gen2, True),
#     "11-2": (rabinowitz_richter.gen3, True),
#     # ERR misprint in Stroud copied from original article
#     # '13-1': (rabinowitz_richter.gen4,
#     "15-1": (rabinowitz_richter.gen5, True),
# }

Stroud = {
    "4-1": stroud_4_1,
    "5-1": StroudSecrest["V"],
    "7-1": StroudSecrest["VI"],
    "9-1": RabinowitzRichter[1],
    "11-1": RabinowitzRichter[2],
    "11-2": RabinowitzRichter[3],
    # ERR misprint in Stroud copied from original article
    # '13-1': RabinowitzRichter[4],
    "15-1": RabinowitzRichter[5],
}
