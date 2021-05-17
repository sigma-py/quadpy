import numpy as np
import sympy

from ..helpers import book
from ._helpers import U2Scheme, register

# Pages 73-74 in
_source = book(
    authors="V.I. Krylov",
    title="Approximate Calculation of Integrals",
    publisher="Macmillan, New York",
    year="1962",
    note="Translated from 1st Russian ed., 1959, by A.H. Stroud",
    url="https://books.google.de/books/about/Approximate_Calculation_of_Integrals.html?id=ELeRwR27IRIC",
)

frac = sympy.Rational
cos = np.vectorize(sympy.cos)
sin = np.vectorize(sympy.sin)
pi = sympy.pi


def krylov(n: int):
    weights = np.full(n, frac(1, n))
    alpha = 2 * np.arange(n) * pi / n
    points = np.column_stack([cos(alpha), sin(alpha)])
    return U2Scheme(f"Krylov {n}", _source, n - 1, weights, points, 2.145e-16)


register([krylov])
