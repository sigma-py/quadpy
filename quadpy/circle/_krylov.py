# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ._helpers import CircleScheme
from ..helpers import book

# Pages 73-74 in
_citation = book(
    authors="V.I. Krylov",
    title="Approximate Calculation of Integrals",
    publisher="Macmillan, New York",
    year="1962",
    note="Translated from 1st Russian ed., 1959, by A.H. Stroud",
    url="https://books.google.de/books/about/Approximate_Calculation_of_Integrals.html?id=ELeRwR27IRIC",
)


def krylov(n, symbolic=False):
    cos = numpy.vectorize(sympy.cos) if symbolic else numpy.cos
    sin = numpy.vectorize(sympy.sin) if symbolic else numpy.sin
    pi = sympy.pi if symbolic else numpy.pi

    weights = numpy.full(n, 2 * pi / n)
    alpha = 2 * numpy.arange(n) * pi / n
    points = numpy.column_stack([cos(alpha), sin(alpha)])
    return CircleScheme("Krylov {}".format(n), _citation, n - 1, weights, points)
