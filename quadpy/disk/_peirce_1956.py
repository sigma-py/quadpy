# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ..helpers import fsd, pm, untangle, fs_array, phdthesis
from ._helpers import DiskScheme


_citation = phdthesis(
    authors=["William Hollis Peirce"],
    title="Numerical integration over planar regions",
    school="University of Wisconsin--Madison",
    year="1956",
    url="https://books.google.de/books/about/Numerical_integration_over_planar_region.html?id=WR9SAAAAMAAJ",
)


def peirce_1956_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pm_ = numpy.array([+1, -1])

    sqrt29 = sqrt(29)
    r = sqrt(frac(3, 4))
    s, t = sqrt((27 - pm_ * 3 * sqrt29) / 104)

    B1 = frac(2, 27)
    # ERR Stroud falsely lists 4 instead of 41 here.
    B2, B3 = (551 + pm_ * 41 * sqrt29) / 6264

    data = [(B1, fsd(2, (r, 1))), (B2, pm(2, s)), (B3, pm(2, t))]

    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Peirce 1956-1", weights, points, 7, _citation)


def peirce_1956_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    cos = numpy.vectorize(sympy.cos) if symbolic else numpy.cos
    sin = numpy.vectorize(sympy.sin) if symbolic else numpy.sin
    pm_ = numpy.array([+1, -1])

    sqrt15 = sqrt(15)
    cos_pi8 = cos(pi / 8)
    sin_pi8 = sin(pi / 8)

    r = sqrt((5 + sqrt15) / 10)
    u1 = sqrt((5 - sqrt15) / 10) * cos_pi8
    v1 = sqrt((5 - sqrt15) / 10) * sin_pi8
    u2 = sqrt(frac(1, 2)) * cos_pi8
    v2 = sqrt(frac(1, 2)) * sin_pi8
    u3, v3 = sqrt((5 + sqrt15 + pm_ * sqrt(185 * sqrt15 - 700)) / 20)

    B1 = (12060 - 1440 * sqrt15) / 254088
    B2 = frac(5, 144)
    B3 = frac(1, 18)
    B4 = (5585 + 1440 * sqrt15) / 508176

    data = [
        (B1, fsd(2, (r, 1))),
        (B2, fs_array([u1, v1])),
        (B3, fs_array([u2, v2])),
        (B4, fs_array([u3, v3])),
    ]

    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Peirce 1956-2", weights, points, 9, _citation)


def peirce_1956_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pm_ = numpy.array([+1, -1])

    sqrt15 = sqrt(15)

    B1 = frac(5, 144)
    B2 = (34 - 5 * sqrt15) / 396
    B3 = (4805 - 620 * sqrt15) / 103824
    C1 = (10 + 5 * sqrt15) / 792
    C2 = (2405 + 620 * sqrt15) / 207648
    D = B1

    r1 = sqrt((5 - sqrt15) / 10)
    r2 = sqrt(frac(1, 2))
    r3 = sqrt((5 + sqrt15) / 10)
    u1, v1 = sqrt((5 + pm_ * sqrt(45 - 10 * sqrt15)) / 20)
    u2, v2 = sqrt((5 + sqrt15 + pm_ * 2 * sqrt(40 * sqrt15 - 150)) / 20)
    t = sqrt((5 - sqrt15) / 20)

    data = [
        (B1, fsd(2, (r1, 1))),
        (B2, fsd(2, (r2, 1))),
        (B3, fsd(2, (r3, 1))),
        (C1, fs_array([u1, v1])),
        (C2, fs_array([u2, v2])),
        (D, pm(2, t)),
    ]

    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Peirce 1956-3", weights, points, 11, _citation)
