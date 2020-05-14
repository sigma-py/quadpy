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
    return CnScheme("Thacher", n, weights, points, 2, _source)
