# ENH in the article, most schemes are given only in single precision. quadpy adds
# symbolic expressions

import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import T2Scheme, register

source = article(
    authors=["P. Hillion"],
    title="Numerical Integration on a Triangle",
    journal="International Journal for Numerical Methods in Engineering",
    volume="11",
    pages="797-815",
    year="1977",
    url="https://doi.org/10.1002/nme.1620110504",
)


def hillion_01():
    d = {"centroid": [[1]]}
    return T2Scheme("Hillion 1", d, 1, source)


def hillion_02():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 2)]]}
    return T2Scheme("Hillion 2", d, 2, source)


def hillion_03():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 6)]]}
    return T2Scheme("Hillion 3", d, 2, source)


def hillion_04():
    a0, a1 = ((3 + i * sqrt(3)) / 8 for i in [+1, -1])
    d = {"s2_static": [[frac(1, 9)], [0]], "swap_ab": [[frac(4, 9)], [a0], [a1]]}
    return T2Scheme("Hillion 4", d, 2, source)


def hillion_05():
    a0, a1 = ((3 + i * sqrt(3)) / 8 for i in [+1, -1])
    d = {
        "s2_static": [[frac(1, 9)], [frac(2, 3)]],
        "swap_ab": [[frac(4, 9)], [frac(2, 3) - a0], [frac(2, 3) - a1]],
    }
    return T2Scheme("Hillion 5", d, 2, source)


def hillion_06():
    lm, mu = ((2 + i * sqrt(2 + i * sqrt(3))) / 6 for i in [+1, -1])
    d = {
        "swap_ab": [
            [frac(1, 4), frac(1, 4)],
            [lm, frac(2, 3) - lm],
            [mu, frac(2, 3) - mu],
        ]
    }
    return T2Scheme("Hillion 6", d, 2, source)


def hillion_07():
    pm = np.array([+1, -1])

    a, b = (6 + sqrt(2) + pm * sqrt(6 * (3 + 2 * sqrt(2)))) / 20
    c, d = (6 - sqrt(2) + pm * sqrt(6 * (3 - 2 * sqrt(2)))) / 20
    w1 = (2 - 3 * (b + c)) / 6 / (a + d - b - c)
    w2 = (2 - 3 * (a + d)) / 6 / (b + c - a - d)

    d = {"swap_ab": [[w1, w2], [a, c], [d, b]]}
    return T2Scheme("Hillion 7", d, 3, source, 2.220e-16)


def hillion_08():
    lambda2, lambda3 = ((32 + i * 2 * sqrt(46)) / 105 for i in [+1, -1])
    w1, w2 = (2 * (3266 + i * 19 * sqrt(46)) / 17664 for i in [+1, -1])
    d = {
        "swap_ab": [[2 * frac(25, 384)], [0], [frac(4, 5)]],
        "s2_static": [[w1, w2], [lambda2, lambda3]],
    }
    return T2Scheme("Hillion 8", d, 3, source)


def hillion_09():
    # ERR the article is missing the minus sign
    d = {"centroid": [[-frac(9, 16)]], "d3_aa": [[frac(25, 48)], [frac(1, 5)]]}
    return T2Scheme("Hillion 9", d, 3, source)


def hillion_10():
    lambda1, lambda2 = ((16 + i * 2 * sqrt(14)) / 25 for i in [+1, -1])
    w1, w2 = (2 * (161 + i * 17 * sqrt(14)) / 2688 for i in [+1, -1])
    d = {
        "swap_ab": [[w2, w1], [lambda1, 0], [0, lambda2]],
        "s2_static": [[frac(25, 48)], [frac(2, 5)]],
    }
    return T2Scheme("Hillion 10", d, 3, source)


register(
    [
        hillion_01,
        hillion_02,
        hillion_03,
        hillion_04,
        hillion_05,
        hillion_06,
        hillion_07,
        hillion_08,
        hillion_09,
        hillion_10,
    ]
)
