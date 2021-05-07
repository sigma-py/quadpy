from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, expand_symmetries
from ._helpers import CnScheme

_source = article(
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
    d = {
        "0": [[frac(8 - 5 * n, 9)]],
        "a0": [[frac(5, 18)], [r]],
        "a": [[frac(1, 9 * 2 ** n)], [1]],
    }
    points, weights = expand_symmetries(d, n)
    return CnScheme("Mustard-Lyness-Blatt", n, weights, points, 5, _source, 6.312e-14)
