# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .albrecht_collatz import AlbrechtCollatz
from .hammer_stroud import HammerStroud
from .hammer_wymore import HammerWymore
from .mustard_lyness_blatt import MustardLynessBlatt
from .sadowsky import Sadowsky
from .sarma_stroud import SarmaStroud
from .stroud1967 import Stroud1967
from .tyler import Tyler

from .helpers import pm_rrr

from ..ncube import Ewing
from ..helpers import untangle


class Stroud(object):
    """
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    """

    def __init__(self, index, symbolic=False):

        d = {
            "C3 3-1": (Tyler, [1]),
            "C3 3-2": (ProductGauss, [symbolic]),
            "C3 3-3": (Ewing, [3, symbolic]),
            "C3 3-4": (MustardLynessBlatt, [1, symbolic]),
            "C3 3-5": (MustardLynessBlatt, [2, symbolic]),
            "C3 3-6": (AlbrechtCollatz, [symbolic]),
            "C3 3-7": (MustardLynessBlatt, [3, symbolic]),
            "C3 5-1": (Stroud1967, [symbolic]),
            "C3 5-2": (HammerStroud, ["2-3", symbolic]),
            "C3 5-3": (Tyler, [2, symbolic]),
            "C3 5-4": (MustardLynessBlatt, [4, symbolic]),
            "C3 5-5": (MustardLynessBlatt, [5, symbolic]),
            "C3 5-6": (MustardLynessBlatt, [6, symbolic]),
            "C3 5-7": (MustardLynessBlatt, [7, symbolic]),
            "C3 5-8": (Sadowsky, [symbolic]),
            "C3 7-1a": (HammerStroud, ["5-3a", symbolic]),
            "C3 7-1b": (HammerStroud, ["5-3b", symbolic]),
            "C3 7-2": (HammerWymore, [symbolic]),
            "C3 7-3": (SarmaStroud, [symbolic]),
        }

        fun, args = d[index]
        scheme = fun(*args)

        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return


class ProductGauss(object):
    def __init__(self, symbolic):
        reference_volume = 8

        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        self.degree = 3
        data = [(frac(1, 8), pm_rrr(sqrt(frac(1, 3))))]
        self.points, self.weights = untangle(data)
        self.weights *= reference_volume
        return
