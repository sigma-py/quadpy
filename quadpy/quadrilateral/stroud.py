# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from __future__ import division

import warnings

import numpy
import sympy

from .albrecht_collatz import AlbrechtCollatz
from .burnside import Burnside
from .irwin import Irwin
from .maxwell import Maxwell
from .meister import Meister
from .miller import Miller
from .phillips import Phillips
from .rabinowitz_richter import RabinowitzRichter
from .tyler import Tyler

from .. import ncube
from .helpers import concat, zero, symm_r0, symm_s, symm_s_t


class ProductTrapezoidal(object):
    def __init__(self, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        reference_volume = 4

        self.degree = 1
        self.weights, self.points = symm_s(frac(1, 4), 1)
        self.weights *= reference_volume
        return


class ProductGauss3(object):
    def __init__(self, symbolic):
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        reference_volume = 4

        self.degree = 3
        # ERR misprint in Stroud: sqrt(1/3) vs 1/3
        self.weights, self.points = symm_s([frac(1, 4), sqrt(frac(1, 3))])
        self.weights *= reference_volume
        return


class ProductGauss5(object):
    def __init__(self, symbolic):
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        reference_volume = 4

        self.degree = 5
        r = sqrt(frac(3, 5))
        self.weights, self.points = concat(
            zero(frac(16, 81)), symm_r0([frac(10, 81), r]), symm_s([frac(25, 324), r])
        )
        self.weights *= reference_volume
        return


class ProductGauss7(object):
    def __init__(self, symbolic):
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        reference_volume = 4

        # TODO fix
        warnings.warn("Formula only has degree 1!")
        self.degree = 1

        pm = numpy.array([+1, -1])

        r, s = sqrt((15 - pm * 2 * sqrt(30)) / 35)

        B1, B2 = (59 + pm * 6 * sqrt(30)) / 864
        B3 = frac(49, 864)

        r = sqrt(frac(3, 5))
        self.weights, self.points = concat(
            symm_s([B1, r], [B2, s]), symm_s_t([B3, r, s])
        )

        self.weights *= reference_volume
        return


Stroud = {
    "C2 1-1": lambda symbolic=False: ProductTrapezoidal(symbolic),
    "C2 1-2": lambda symbolic=False: Miller(symbolic),
    "C2 3-1": lambda symbolic=False: ProductGauss3(symbolic),
    "C2 3-2": lambda symbolic=False: ncube.Ewing(2, symbolic),
    # product Simpson:
    "C2 3-3": lambda symbolic=False: ncube.Stroud(2, "Cn 3-6", symbolic),
    "C2 3-4": lambda symbolic=False: AlbrechtCollatz[1](symbolic),
    "C2 3-5": lambda symbolic=False: Irwin[1](symbolic),
    "C2 5-1": lambda symbolic=False: AlbrechtCollatz[2](symbolic),
    "C2 5-2": lambda symbolic=False: AlbrechtCollatz[3](symbolic),
    "C2 5-3": lambda symbolic=False: Burnside(symbolic),
    "C2 5-4": lambda symbolic=False: ProductGauss5(symbolic),
    "C2 5-5": lambda symbolic=False: Tyler(1, symbolic),
    "C2 5-6": lambda symbolic=False: AlbrechtCollatz[4](symbolic),
    "C2 5-7": lambda symbolic=False: Irwin[2](symbolic),
    "C2 7-1": lambda symbolic=False: Tyler(2, symbolic),
    "C2 7-2": lambda symbolic=False: Phillips(symbolic),
    "C2 7-3": lambda symbolic=False: Maxwell(symbolic),
    "C2 7-4": lambda symbolic=False: ProductGauss7(symbolic),
    "C2 7-5": lambda symbolic=False: Tyler(3, symbolic),
    "C2 7-6": lambda symbolic=False: Meister(symbolic),
    "C2 9-1": RabinowitzRichter[1],
    "C2 11-1": RabinowitzRichter[2],
    "C2 11-2": RabinowitzRichter[3],
    "C2 13-1": RabinowitzRichter[4],
    "C2 15-1": RabinowitzRichter[5],
    "C2 15-2": RabinowitzRichter[6],
}

# elif index == 2:
#     self.weights = 4 * [1.0]
#     self.points = _symm_r_0(sqrt(2.0/3.0))
#     self.degree = 3
# elif index == 3:
#     self.weights = (
#         4 * [2.0/3.0] +
#         [4.0/3.0]
#         )
#     self.points = (
#         _symm_r_0(1.0) +
#         [[0.0, 0.0]]
#         )
#     self.degree = 3
# elif index == 4:
#     self.weights = (
#         4 * [1.0/3.0] +
#         [8.0/3.0]
#         )
#     self.points = (
#         _symm_s(1.0) +
#         [[0.0, 0.0]]
#         )
#     self.degree = 3
# elif index == 6:
#     self.weights = (
#         4 * [1.0 / 9.0] +
#         1 * [-8.0 / 9.0] +
#         4 * [10.0 / 9.0]
#         )
#     self.points = (
#         _symm_s(1.0) +
#         [[0.0, 0.0]] +
#         _symm_r_0(sqrt(2.0 / 5.0))
#         )
#     self.degree = 5
# elif index == 7:
#     self.weights = (
#         [8.0 / 7.0] +
#         8 * [5.0 / 14.0]
#         )
#     self.points = (
#         [[0.0, 0.0]] +
#         _symm_s_t(
#             0.846233119448533574334773553209421,
#             0.466071712166418914664132375714293
#             )
#         )
#     self.degree = 5
# else:
#     assert index == 10
#     scheme1d = line_segment.GaussLegendre(8)
#     self.weights = numpy.outer(
#         scheme1d.weights, scheme1d.weights
#         ).flatten()
#     self.points = numpy.dstack(numpy.meshgrid(
#         scheme1d.points, scheme1d.points
#         )).reshape(-1, 2)
#     assert len(self.points) == 64
#     self.degree = 15
