# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .hammer_stroud import HammerStroud
from .lauffer import Lauffer
from .stroud1961 import Stroud1961
from .stroud1964 import Stroud1964
from .stroud1966 import Stroud1966
from .stroud1969 import Stroud1969

from ..helpers import untangle


class Stroud(object):
    """
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    """

    def __init__(self, n, index, symbolic=False):
        self.name = "Stroud({})".format(index)
        self.dim = n

        scheme = {
            "Tn 1-1": lambda: MidPoint(n, symbolic),
            "Tn 1-2": lambda: Lauffer(n, 1, symbolic),
            "Tn 2-1a": lambda: HammerStroud(n, "1a", symbolic),
            "Tn 2-1b": lambda: HammerStroud(n, "1b", symbolic),
            "Tn 2-2": lambda: Lauffer(n, 2, symbolic),
            "Tn 3-1": lambda: HammerStroud(n, "2", symbolic),
            "Tn 3-2": lambda: Stroud1966(n, "I", symbolic),
            "Tn 3-3": lambda: Stroud1961(n, symbolic),
            "Tn 3-4": lambda: Stroud1966(n, "II", symbolic),
            "Tn 3-5": lambda: Stroud1966(n, "III", symbolic),
            "Tn 3-6a": lambda: Stroud1964(n, "a", symbolic),
            "Tn 3-6b": lambda: Stroud1964(n, "b", symbolic),
            "Tn 3-7": lambda: Stroud1966(n, "IV", symbolic),
            "Tn 3-8": lambda: Stroud1966(n, "V", symbolic),
            "Tn 3-9": lambda: Lauffer(n, 3, symbolic),
            "Tn 3-10": lambda: Stroud1966(n, "VI", symbolic),
            "Tn 3-11": lambda: Stroud1966(n, "VII", symbolic),
            "Tn 4-1": lambda: Lauffer(n, 4, symbolic),
            "Tn 5-1": lambda: Stroud1969(n, symbolic),
            "Tn 5-2": lambda: Lauffer(n, 5, symbolic),
        }[index]()

        self.degree = scheme.degree
        self.weights = scheme.weights
        self.bary = scheme.bary
        self.points = scheme.points
        return


class MidPoint(object):
    def __init__(self, n, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.degree = 1
        data = [(1, numpy.full((1, n + 1), frac(1, n + 1)))]
        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
