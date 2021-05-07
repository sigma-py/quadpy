import numpy as np
from sympy import sqrt

from ..helpers import article, untangle
from ._helpers import CnScheme, _s

_source = article(
    authors=["Henry C. Thacher"],
    title="An efficient composite formula for multidimensional quadrature",
    journal="Communications of the ACM",
    volume="7",
    number="1",
    month="jan",
    year="1964",
    pages="23-25",
    url="https://doi.org/10.1145/363872.363897",
)


def thacher(n):
    r = sqrt(3) / 6
    data = [(1, [n * [2 * r]]), (+r, _s(n, -1, r)), (-r, _s(n, +1, r))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return CnScheme(f"Thacher ({n}D)", n, weights, points, 2, _source, 1.399e-14)
