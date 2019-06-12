# -*- coding: utf-8 -*-
#
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
from .helpers import _symm_r_0, _symm_s, _symm_s_t, _z

from .. import ncube
from ..helpers import untangle


class Stroud(object):
    """
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    """

    def __init__(self, index, symbolic=False):
        self.name = "Stroud({})".format(index)

        scheme = {
            "C2 1-1": lambda: ProductTrapezoidal(symbolic),
            "C2 1-2": lambda: Miller(symbolic),
            "C2 3-1": lambda: ProductGauss3(symbolic),
            "C2 3-2": lambda: ncube.Ewing(2, symbolic),
            # product Simpson:
            "C2 3-3": lambda: ncube.Stroud(2, "Cn 3-6", symbolic),
            "C2 3-4": lambda: AlbrechtCollatz[1](symbolic),
            "C2 3-5": lambda: Irwin(1, symbolic),
            "C2 5-1": lambda: AlbrechtCollatz[2](symbolic),
            "C2 5-2": lambda: AlbrechtCollatz[3](symbolic),
            "C2 5-3": lambda: Burnside(symbolic),
            "C2 5-4": lambda: ProductGauss5(symbolic),
            "C2 5-5": lambda: Tyler(1, symbolic),
            "C2 5-6": lambda: AlbrechtCollatz[4](symbolic),
            "C2 5-7": lambda: Irwin(2, symbolic),
            "C2 7-1": lambda: Tyler(2, symbolic),
            "C2 7-2": lambda: Phillips(symbolic),
            "C2 7-3": lambda: Maxwell(symbolic),
            "C2 7-4": lambda: ProductGauss7(symbolic),
            "C2 7-5": lambda: Tyler(3, symbolic),
            "C2 7-6": lambda: Meister(symbolic),
            "C2 9-1": lambda: RabinowitzRichter(1),
            "C2 11-1": lambda: RabinowitzRichter(2),
            "C2 11-2": lambda: RabinowitzRichter(3),
            "C2 13-1": lambda: RabinowitzRichter(4),
            "C2 15-1": lambda: RabinowitzRichter(5),
            "C2 15-2": lambda: RabinowitzRichter(6),
        }[index]()

        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points

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
        return


class ProductTrapezoidal(object):
    def __init__(self, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        reference_volume = 4

        self.degree = 1
        self.weights = reference_volume * numpy.full(4, frac(1, 4))
        self.points = _symm_s(1)
        return


class ProductGauss3(object):
    def __init__(self, symbolic):
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        reference_volume = 4

        self.degree = 3
        self.weights = reference_volume * numpy.full(4, frac(1, 4))
        # ERR misprint in Stroud: sqrt(1/3) vs 1/3
        self.points = _symm_s(sqrt(frac(1, 3)))
        return


class ProductGauss5(object):
    def __init__(self, symbolic):
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        reference_volume = 4

        self.degree = 5
        r = sqrt(frac(3, 5))
        data = [
            (frac(16, 81), _z()),
            (frac(10, 81), _symm_r_0(r)),
            (frac(25, 324), _symm_s(r)),
        ]
        self.points, self.weights = untangle(data)
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
        data = [(B1, _symm_s(r)), (B2, _symm_s(s)), (B3, _symm_s_t(r, s))]
        self.points, self.weights = untangle(data)
        self.weights *= reference_volume
        return
