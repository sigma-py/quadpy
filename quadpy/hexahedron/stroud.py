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

from .albrecht_collatz import albrecht_collatz
from .hammer_stroud import hammer_stroud_2_3, hammer_stroud_5_3a, hammer_stroud_5_3b
from .hammer_wymore import hammer_wymore
from .mustard_lyness_blatt import (
    mustard_lyness_blatt_1,
    mustard_lyness_blatt_2,
    mustard_lyness_blatt_3,
    mustard_lyness_blatt_4,
    mustard_lyness_blatt_5,
    mustard_lyness_blatt_6,
    mustard_lyness_blatt_7,
)
from .sadowsky import sadowsky
from .sarma_stroud import sarma_stroud
from .stroud1967 import stroud_1967
from .tyler import tyler_1, tyler_2

from .helpers import pm_rrr

from ..ncube import Ewing
from ..helpers import untangle
from .helpers import HexahedronScheme


def stroud_3_2(symbolic=False):
    # Product Gauss scheme
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    data = [(frac(1, 8), pm_rrr(sqrt(frac(1, 3))))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Stroud 3-2", 3, weights, points)


Stroud = {
    "C3 3-1": tyler_1,
    "C3 3-2": stroud_3_2,
    "C3 3-3": lambda symbolic=False: Ewing(3, symbolic),
    "C3 3-4": mustard_lyness_blatt_1,
    "C3 3-5": mustard_lyness_blatt_2,
    "C3 3-6": albrecht_collatz,
    "C3 3-7": mustard_lyness_blatt_3,
    "C3 5-1": stroud_1967,
    "C3 5-2": hammer_stroud_2_3,
    "C3 5-3": tyler_2,
    "C3 5-4": mustard_lyness_blatt_4,
    "C3 5-5": mustard_lyness_blatt_5,
    "C3 5-6": mustard_lyness_blatt_6,
    "C3 5-7": mustard_lyness_blatt_7,
    "C3 5-8": sadowsky,
    "C3 7-1a": hammer_stroud_5_3a,
    "C3 7-1b": hammer_stroud_5_3b,
    "C3 7-2": hammer_wymore,
    "C3 7-3": sarma_stroud,
}
