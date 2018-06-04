# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .ewing import Ewing
from .hammer_stroud import HammerStroud
from .mustard_lyness_blatt import MustardLynessBlatt
from .phillips import Phillips
from .stroud1957 import Stroud1957
from .stroud1966 import Stroud1966
from .stroud1968 import Stroud1968
from .thacher import Thacher
from .tyler import Tyler

from ..helpers import fsd, pm


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
            "Cn 1-1": lambda: Centroid(n),
            "Cn 1-2": lambda: ProductTrapezoidal(n),
            "Cn 2-1": lambda: Stroud1957(n, 2, symbolic),
            "Cn 2-2": lambda: Thacher(n, symbolic),
            "Cn 3-1": lambda: Stroud1957(n, 3, symbolic),
            "Cn 3-2": lambda: Cn32(n, symbolic),
            "Cn 3-3": lambda: Tyler(n, symbolic),
            "Cn 3-4": lambda: ProductGauss3(n, symbolic),
            "Cn 3-5": lambda: Ewing(n, symbolic),
            "Cn 3-6": lambda: ProductSimpson(n, symbolic),
            # TODO implement Cn 5-1
            # Cn 5-1 is not implemented because it's based on explicit values
            # only given for n=4,5,6.
            "Cn 5-2": lambda: HammerStroud(n, "2-n", symbolic),
            "Cn 5-3": lambda: Stroud1968(n, symbolic),
            "Cn 5-4": lambda: Stroud1966(n, "a", symbolic),
            "Cn 5-5": lambda: MustardLynessBlatt(n, symbolic),
            "Cn 5-6": lambda: Stroud1966(n, "b", symbolic),
            "Cn 5-7": lambda: Stroud1966(n, "c", symbolic),
            "Cn 5-8": lambda: Stroud1966(n, "d", symbolic),
            "Cn 5-9": lambda: ProductGauss5(n, symbolic),
            "Cn 7-1": lambda: Phillips(n, symbolic),
        }[index]()

        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return


class Centroid(object):
    def __init__(self, n):
        reference_volume = 2 ** n
        self.degree = 1
        self.weights = numpy.array([reference_volume])
        self.points = numpy.full((1, n), 0)
        return


class ProductTrapezoidal(object):
    def __init__(self, n):
        self.degree = 1
        self.weights = numpy.full(2 ** n, 1)
        self.points = pm(n, 1)
        return


class Cn32(object):
    def __init__(self, n, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        reference_volume = 2 ** n

        self.degree = 3
        self.weights = numpy.full(2 * n, frac(reference_volume, 2 * n))
        r = sqrt(frac(n, 3))
        self.points = fsd(n, (r, 1))
        return


class ProductGauss3(object):
    def __init__(self, n, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 3
        reference_volume = 2 ** n
        self.weights = numpy.full(2 ** n, frac(reference_volume, 2 ** n))
        r = sqrt(3) / 3
        self.points = pm(n, r)
        return


class ProductGauss5(object):
    def __init__(self, n, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 5
        lst = n * [[frac(5, 9), frac(8, 9), frac(5, 9)]]
        self.weights = numpy.product(
            numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n), axis=-1
        )
        sqrt35 = sqrt(frac(3, 5))
        lst = n * [[-sqrt35, 0, sqrt35]]
        self.points = numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n)
        return


class ProductSimpson(object):
    def __init__(self, n, symbolic):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.degree = 3
        lst = n * [[frac(1, 3), frac(4, 3), frac(1, 3)]]
        self.weights = numpy.product(
            numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n), axis=-1
        )
        lst = n * [[-1, 0, +1]]
        self.points = numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n)
        return
