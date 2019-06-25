# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ._helpers import TriangleScheme
from ..helpers import techreport

citation = techreport(
    authors=["Noel J. Walkington"],
    title="Quadrature on simplices of arbitrary dimension",
    institution="CMU",
    year="2000",
    url="https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/",
)


def walkington_p5(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 5

    a1, a2 = [(155 + i * sqrt(15)) / 1200 for i in [+1, -1]]
    weights = numpy.concatenate(
        [numpy.full(1, frac(9, 40)), numpy.full(3, a1), numpy.full(3, a2)]
    )

    x1, x2 = [(6 + i * sqrt(15)) / 21 for i in [+1, -1]]
    points = numpy.concatenate([_c(frac), _xi1(x1), _xi1(x2)])
    return TriangleScheme("Walkington p5", weights, points, degree, citation)


def _c(frac):
    return numpy.full((1, 3), frac(1, 3))


def _xi1(a):
    b = 1 - 2 * a
    return numpy.array([[b, a, a], [a, b, a], [a, a, b]])
