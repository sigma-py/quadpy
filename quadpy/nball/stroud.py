# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .hammer_stroud import HammerStroud
from .stenger import Stenger
from .stroud1957 import Stroud1957
from .stroud1966 import Stroud1966
from .stroud1967a import Stroud1967a
from .stroud1967b import Stroud1967b

from .helpers import volume_unit_ball
from ..helpers import pm, untangle


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
            "Sn 2-1": lambda: Stroud1957(n, symbolic),
            "Sn 3-1": lambda: HammerStroud(n, "11-n", 0, symbolic),
            "Sn 3-2": lambda: Sn32(n, symbolic),
            "Sn 5-1a": lambda: Stroud1967a(n, "a"),
            "Sn 5-1b": lambda: Stroud1967a(n, "b"),
            "Sn 5-2": lambda: HammerStroud(n, "12-n", 0, symbolic),
            "Sn 5-3": lambda: Stroud1966(n, "a", symbolic),
            "Sn 5-4": lambda: Stroud1966(n, "b", symbolic),
            "Sn 5-5": lambda: Stroud1966(n, "c", symbolic),
            "Sn 5-6": lambda: Stroud1966(n, "d", symbolic),
            "Sn 7-1a": lambda: Stroud1967b(n, "a", symbolic),
            "Sn 7-1b": lambda: Stroud1967b(n, "b", symbolic),
            "Sn 7-2": lambda: Stroud1967b(n, "c", symbolic),
            # TODO fix
            # "Sn 7-3a": lambda: Stenger(n, 7, "a"),
            "Sn 7-3b": lambda: Stenger(n, 7, "b"),
            "Sn 9-1a": lambda: Stenger(n, 9, "a"),
            "Sn 9-1b": lambda: Stenger(n, 9, "b"),
            "Sn 11-1a": lambda: Stenger(n, 11, "a"),
            "Sn 11-1b": lambda: Stenger(n, 11, "b"),
        }[index]()

        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return


class Sn32(object):
    def __init__(self, n, symbolic):
        self.degree = 3

        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        r = sqrt(frac(1, n + 2))
        data = [(frac(1, 2 ** n), pm(n, r))]
        self.points, self.weights = untangle(data)
        self.weights *= volume_unit_ball(n, symbolic=symbolic)
        return
