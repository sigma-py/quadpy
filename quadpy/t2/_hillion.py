# ENH in the article, most schemes are given only in single precision. quadpy adds
# symbolic expressions

import numpy
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import T2Scheme, expand_symmetries

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
    d = {"s3": [[frac(1, 2)]]}
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 1", weights, points, 1, source)


def hillion_02():
    d = {"s2": [[frac(1, 6)], [frac(1, 2)]]}
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 2", weights, points, 2, source)


def hillion_03():
    d = {"s2": [[frac(1, 6)], [frac(1, 6)]]}
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 3", weights, points, 2, source)


def hillion_04():
    a0, a1 = [(3 + i * sqrt(3)) / 8 for i in [+1, -1]]
    d = {"s2_static": [[frac(1, 18)], [0]], "swap_ab": [[frac(2, 9)], [a0], [a1]]}
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 4", weights, points, 2, source)


def hillion_05():
    a0, a1 = [(3 + i * sqrt(3)) / 8 for i in [+1, -1]]
    d = {
        "s2_static": [[frac(1, 18)], [frac(2, 3)]],
        "swap_ab": [[frac(2, 9)], [frac(2, 3) - a0], [frac(2, 3) - a1]],
    }
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 5", weights, points, 2, source)


def hillion_06():
    lm, mu = [(2 + i * sqrt(2 + i * sqrt(3))) / 6 for i in [+1, -1]]
    d = {
        "swap_ab": [
            [frac(1, 8), frac(1, 8)],
            [lm, frac(2, 3) - lm],
            [mu, frac(2, 3) - mu],
        ]
    }
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 6", weights, points, 2, source)


def hillion_07():
    pm = numpy.array([+1, -1])

    a, b = (6 + sqrt(2) + pm * sqrt(6 * (3 + 2 * sqrt(2)))) / 20
    c, d = (6 - sqrt(2) + pm * sqrt(6 * (3 - 2 * sqrt(2)))) / 20
    w1 = (2 - 3 * (b + c)) / 12 / (a + d - b - c)
    w2 = (2 - 3 * (a + d)) / 12 / (b + c - a - d)

    d = {"swap_ab": [[w1, w2], [a, c], [d, b]]}
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 7", weights, points, 3, source)


def hillion_08():
    lambda2, lambda3 = [(32 + i * 2 * sqrt(46)) / 105 for i in [+1, -1]]
    w1, w2 = [(3266 + i * 19 * sqrt(46)) / 17664 for i in [+1, -1]]
    d = {
        "swap_ab": [[frac(25, 384)], [0], [frac(4, 5)]],
        "s2_static": [[w1, w2], [lambda2, lambda3]],
    }
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 8", weights, points, 3, source)


def hillion_09():
    # ERR the article is missing the minus sign
    d = {"s3": [[-frac(9, 32)]], "s2": [[frac(25, 96)], [frac(1, 5)]]}
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 9", weights, points, 3, source)


def hillion_10():
    lambda1, lambda2 = [(16 + i * 2 * sqrt(14)) / 25 for i in [+1, -1]]
    w1, w2 = [(161 + i * 17 * sqrt(14)) / 2688 for i in [+1, -1]]
    d = {
        "swap_ab": [[w2, w1], [lambda1, 0], [0, lambda2]],
        "s2_static": [[frac(25, 96)], [frac(2, 5)]],
    }
    points, weights = expand_symmetries(d)
    weights *= 2
    return T2Scheme("Hillion 10", weights, points, 3, source)
