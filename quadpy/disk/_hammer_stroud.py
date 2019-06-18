# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy
import numpy

from ..helpers import untangle, fsd, pm, z, fs_array as fs, article
from .helpers import DiskScheme

from ._radon import radon
from ._peirce_1956 import peirce_1956_1, peirce_1956_3


_citation = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    volume="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)


def hammer_stroud_11_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    # ERR Wrongly stated in Stroud with 0.5 instead of sqrt(0.5)
    data = [(frac(1, 4), fsd(2, (sqrt(frac(1, 2)), 1)))]
    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Hammer-Stroud 11-2", weights, points, 3, _citation)


def hammer_stroud_12_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    data = [
        (frac(1, 6), z(2)),
        (frac(1, 6), fsd(2, (sqrt(frac(1, 2)), 1))),
        (frac(1, 24), pm(2, sqrt(frac(1, 2)))),
    ]
    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Hammer-Stroud 12-2", weights, points, 5, _citation)


def hammer_stroud_13_2(symbolic=False):
    return peirce_1956_1(symbolic)


def hammer_stroud_17(symbolic=False):
    # ENH This is Radon's formula.
    return radon(0, symbolic)


def hammer_stroud_18(symbolic=False):
    # ENH The article only gives floats, but really this is the spherical-product gauss
    # formula as described in Strouds book, S2 7-2.
    #
    # data = [
    #     (frac(1, 16), fs([0.4247082002778669, 0.1759198966061612])),
    #     (frac(1, 16), fs([0.8204732385702833, 0.3398511429799874])),
    # ]
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pm_ = numpy.array([+1, -1])
    cos = numpy.vectorize(sympy.cos) if symbolic else numpy.cos
    sin = numpy.vectorize(sympy.sin) if symbolic else numpy.sin
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi

    r1, r2 = sqrt((3 - pm_ * sqrt(3)) / 6)

    a = (2 * numpy.arange(8) + 1) * pi / 8
    x = numpy.array([cos(a), sin(a)]).T

    data = [(frac(1, 16), r1 * x), (frac(1, 16), r2 * x)]
    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Hammer-Stroud 18", weights, points, 7, _citation)


def hammer_stroud_19(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt6 = sqrt(6)
    alpha1 = (16 + sqrt6) / 288
    alpha2 = (137 - 32 * sqrt6) / 1818
    alpha3 = (520 + 155 * sqrt6) / 3636 / 8

    data = [
        (frac(1, 9), z(2)),
        (alpha1, fs([0.5505043204538557, 0.2280263556769715])),
        (alpha2, fsd(2, (0.9192110607898046, 1))),
        (alpha3, fs([0.7932084745126058, 0.4645097310495256])),
    ]
    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Hammer-Stroud 19", weights, points, 9, _citation)


def hammer_stroud_20(symbolic=False):
    # ENH Also Peirce's formula, even given symbolically.
    return peirce_1956_3(symbolic)


def hammer_stroud_21(symbolic=False):
    pi = sympy.pi if symbolic else numpy.pi

    alpha0 = 0.0341505695624825 / pi
    alpha1 = 0.0640242008621985 / pi
    alpha2 = 0.0341505695624825 / pi

    data = [
        (alpha0, fs([0.2584361661674054, 0.0514061496288813])),
        (alpha1, fs([0.5634263397544869, 0.1120724670846205])),
        (alpha1, fs([0.4776497869993547, 0.3191553840796721])),
        (alpha1, fs([0.8028016728473508, 0.1596871812824163])),
        (alpha1, fs([0.6805823955716280, 0.4547506180649039])),
        (alpha2, fs([0.2190916025980981, 0.1463923286035535])),
        (alpha2, fs([0.9461239423417719, 0.1881957532057769])),
        (alpha2, fs([0.8020851487551318, 0.5359361621905023])),
    ]

    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Hammer-Stroud 21", weights, points, 15, _citation)
