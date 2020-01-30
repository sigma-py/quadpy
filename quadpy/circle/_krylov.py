import numpy
import sympy

from ..helpers import book
from ._helpers import CircleScheme

# Pages 73-74 in
_citation = book(
    authors="V.I. Krylov",
    title="Approximate Calculation of Integrals",
    publisher="Macmillan, New York",
    year="1962",
    note="Translated from 1st Russian ed., 1959, by A.H. Stroud",
    url="https://books.google.de/books/about/Approximate_Calculation_of_Integrals.html?id=ELeRwR27IRIC",
)

cos = numpy.vectorize(sympy.cos)
sin = numpy.vectorize(sympy.sin)
pi = sympy.pi


def krylov(n):
    weights = numpy.full(n, 2 * pi / n)
    alpha = 2 * numpy.arange(n) * pi / n
    points = numpy.column_stack([cos(alpha), sin(alpha)])
    return CircleScheme(f"Krylov {n}", _citation, n - 1, weights, points)
