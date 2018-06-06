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

        # Adding lambdas here makes sure that only the selected scheme is actually
        # created.
        scheme = {
            "C3 3-1": lambda: Tyler(1),
            "C3 3-2": lambda: ProductGauss(symbolic=symbolic),
            "C3 3-3": lambda: Ewing(3, symbolic=symbolic),
            "C3 3-4": lambda: MustardLynessBlatt(1, symbolic=symbolic),
            "C3 3-5": lambda: MustardLynessBlatt(2, symbolic=symbolic),
            "C3 3-6": lambda: AlbrechtCollatz(symbolic=symbolic),
            "C3 3-7": lambda: MustardLynessBlatt(3, symbolic=symbolic),
            "C3 5-1": lambda: Stroud1967(symbolic=symbolic),
            "C3 5-2": lambda: HammerStroud("2-3", symbolic=symbolic),
            "C3 5-3": lambda: Tyler(2, symbolic=symbolic),
            "C3 5-4": lambda: MustardLynessBlatt(4, symbolic=symbolic),
            "C3 5-5": lambda: MustardLynessBlatt(5, symbolic=symbolic),
            "C3 5-6": lambda: MustardLynessBlatt(6, symbolic=symbolic),
            "C3 5-7": lambda: MustardLynessBlatt(7, symbolic=symbolic),
            "C3 5-8": lambda: Sadowsky(symbolic=symbolic),
            "C3 7-1a": lambda: HammerStroud("5-3a", symbolic=symbolic),
            "C3 7-1b": lambda: HammerStroud("5-3b", symbolic=symbolic),
            "C3 7-2": lambda: HammerWymore(symbolic=symbolic),
            "C3 7-3": lambda: SarmaStroud(symbolic=symbolic),
        }[index]()

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
