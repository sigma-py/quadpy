from math import factorial

import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import techreport, untangle
from ._helpers import TnScheme

source = techreport(
    authors=["Noel J. Walkington"],
    title="Quadrature on simplices of arbitrary dimension",
    institution="CMU",
    year="2000",
    url="https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/",
)


def walkington_1(d):
    degree = 1
    data = [(frac(1, factorial(d)), _c(d, frac))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    # normalize weights
    weights /= np.sum(weights)
    return TnScheme("Walkington 1", d, weights, points, degree, source)


def walkington_2(d):
    # ERR The article claims order 2, but tests really only show order 1.
    # Also, the article says:
    #
    # > The points are inside the simplex when the positive square root is selected.
    #
    # Not sure what this means, but for d>=2, the points are outside the simplex.
    degree = 1
    data = [(frac(1, factorial(d + 1)), _xi1(d, 1 / sqrt(d + 1)))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    # normalize weights
    weights /= np.sum(weights)
    return TnScheme("Walkington 2", d, weights, points, degree, source)


def walkington_3(d):
    degree = 3
    data = [
        (frac(-((d + 1) ** 3), 4 * factorial(d + 2)), _c(d, frac)),
        (frac(+((d + 3) ** 3), 4 * factorial(d + 3)), _xi1(d, frac(1, (d + 3)))),
    ]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    # normalize weights
    weights /= np.sum(weights)
    return TnScheme("Walkington 3", d, weights, points, degree, source)


def walkington_5(d):
    degree = 5
    w0 = frac(+((d + 1) ** 5), 32 * factorial(d + 3))
    w1 = frac(-((d + 3) ** 5), 16 * factorial(d + 4))
    w2 = frac(+((d + 5) ** 5), 16 * factorial(d + 5))
    data = [
        (w0, _c(d, frac)),
        (w1, _xi1(d, frac(1, d + 3))),
        (w2, _xi1(d, frac(1, d + 5))),
        (w2, _xi11(d, frac(1, d + 5), frac)),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    # normalize weights
    weights /= np.sum(weights)
    return TnScheme("Walkington 5", d, weights, points, degree, source)


def walkington_7(d):
    degree = 7
    w0 = -frac(1, 384) * frac((d + 1) ** 7, factorial(d + 4))
    w1 = +frac(1, 128) * frac((d + 3) ** 7, factorial(d + 5))
    w2 = -frac(1, 64) * frac((d + 5) ** 7, factorial(d + 6))
    w3 = +frac(1, 64) * frac((d + 7) ** 7, factorial(d + 7))
    data = [
        (w0, _c(d, frac)),
        (w1, _xi1(d, frac(1, d + 3))),
        (w2, _xi1(d, frac(1, d + 5))),
        (w2, _xi11(d, frac(1, d + 5), frac)),
        (w3, _xi1(d, frac(1, d + 7))),
        (w3, _xi21(d, frac(1, d + 7), frac)),
        (w3, _xi111(d, frac(1, d + 7), frac)),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    # normalize weights
    weights /= np.sum(weights)
    return TnScheme("Walkington 7", d, weights, points, degree, source)


def _c(d, frac):
    return np.array([np.full(d + 1, frac(1, d + 1))])


def _xi1(d, a):
    out = np.full((d + 1, d + 1), a)
    b = 1 - d * a
    np.fill_diagonal(out, b)
    return out


def _xi11(d, a, frac):
    assert d > 1
    b = frac(1 - (d - 1) * a, 2)
    if d == 2:
        out = np.array([[b, b, a], [b, a, b], [a, b, b]])
    else:
        assert d == 3
        out = np.array(
            [
                [b, b, a, a],
                [b, a, b, a],
                [b, a, a, b],
                [a, b, a, b],
                [a, a, b, b],
                [a, b, b, a],
            ]
        )
    return out


def _xi21(d, a, frac):
    assert d > 1
    b = frac(1 - (d - 2) * a, 3)
    # ERR Note that the article incorrectly states (d-2) the the expression for c.
    c = 1 - (d - 1) * a - b
    if d == 2:
        out = np.array(
            [[b, c, a], [c, b, a], [c, a, b], [b, a, c], [a, b, c], [a, c, b]]
        )
    else:
        assert d == 3
        out = np.array(
            [
                [b, c, a, a],
                [b, a, c, a],
                [b, a, a, c],
                [a, b, a, c],
                [a, a, b, c],
                [a, b, c, a],
                [c, b, a, a],
                [c, a, b, a],
                [c, a, a, b],
                [a, c, a, b],
                [a, a, c, b],
                [a, c, b, a],
            ]
        )

    return out


def _xi111(d, a, frac):
    assert d == 3
    b = frac(1 - (d - 2) * a, 3)
    out = np.array([[b, b, b, a], [b, b, a, b], [b, a, b, b], [a, b, b, b]])
    return out
