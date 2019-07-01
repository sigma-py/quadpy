from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, pm, untangle, z
from ._helpers import NCubeScheme

_citation = article(
    authors=["D. Mustard", "J.N. Lyness", "J.M. Blatt"],
    title="Numerical quadrature in n dimensions",
    journal="Comput J",
    year="1963",
    volume="6",
    number="1",
    pages="75-87",
    url="https://doi.org/10.1093/comjnl/6.1.75",
)


def mustard_lyness_blatt(n):
    r = sqrt(frac(2, 5))
    data = [
        (frac(8 - 5 * n, 9), z(n)),
        (frac(5, 18), fsd(n, (r, 1))),
        (frac(1, 9 * 2 ** n), pm(n, 1)),
    ]

    points, weights = untangle(data)
    weights *= 2 ** n
    return NCubeScheme("Mustard-Lyness-Blatt", n, weights, points, 5, _citation)
