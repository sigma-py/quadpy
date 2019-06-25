# -*- coding: utf-8 -*-
#
import numpy
import sympy

from ._helpers import BallScheme
from ..helpers import untangle, fsd, pm, article

_citation = article(
    authors=["I.P. Mysovskih"],
    title="On the construction of cubature formulas for the simplest regions",
    journal="Z. Vychisl. Mat. i. Mat. Fiz.",
    number="4",
    pages="3-14",
    year="1964",
)


def mysovskih(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi
    pm_ = numpy.array([+1, -1])

    sqrt17770 = sqrt(17770)
    r, s = sqrt((1715 - pm_ * 7 * sqrt17770) / 2817)
    t = sqrt(frac(7, 18))
    u = sqrt(frac(7, 27))

    B1, B2 = (2965 * sqrt17770 + pm_ * 227816) / 72030 / sqrt17770
    B3 = frac(324, 12005)
    B4 = frac(2187, 96040)

    data = [
        (B1, fsd(3, (r, 1))),
        (B2, fsd(3, (s, 1))),
        (B3, fsd(3, (t, 2))),
        (B4, pm(3, u)),
    ]

    points, weights = untangle(data)
    weights *= frac(4, 3) * pi
    return BallScheme("Mysovskih", _citation, 7, weights, points)
