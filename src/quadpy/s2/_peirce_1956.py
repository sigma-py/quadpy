import numpy as np
import sympy

from ..helpers import phdthesis
from ._helpers import S2Scheme, register

_source = phdthesis(
    authors=["William Hollis Peirce"],
    title="Numerical integration over planar regions",
    school="University of Wisconsin--Madison",
    year="1956",
    url="https://books.google.de/books/about/Numerical_integration_over_planar_region.html?id=WR9SAAAAMAAJ",
)

frac = sympy.Rational
pi = sympy.pi
sqrt = np.vectorize(sympy.sqrt)
cos = np.vectorize(sympy.cos)
sin = np.vectorize(sympy.sin)
pm_ = np.array([+1, -1])


def peirce_1956_1():
    sqrt29 = sqrt(29)
    r = sqrt(frac(3, 4))
    s, t = sqrt((27 - pm_ * 3 * sqrt29) / 104)

    B1 = frac(2, 27)
    # ERR Stroud incorrectly lists 4 instead of 41 here.
    B2, B3 = (551 + pm_ * 41 * sqrt29) / 6264

    d = {"d4_a0": [[B1], [r]], "d4_aa": [[B2, B3], [s, t]]}
    return S2Scheme("Peirce 1956-1", d, 7, _source)


def peirce_1956_2():
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

    d = {
        "d4_a0": [[B1], [r]],
        "d4_ab": [[B2, B3, B4], [u1, u2, u3], [v1, v2, v3]],
    }
    return S2Scheme("Peirce 1956-2", d, 9, _source)


def peirce_1956_3():
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

    d = {
        "d4_a0": [[B1, B2, B3], [r1, r2, r3]],
        "d4_ab": [[C1, C2], [u1, u2], [v1, v2]],
        "d4_aa": [[D], [t]],
    }
    return S2Scheme("Peirce 1956-3", d, 11, _source)


register([peirce_1956_1, peirce_1956_2, peirce_1956_3])
