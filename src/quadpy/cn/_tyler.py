from sympy import Rational as frac

from ..helpers import article, expand_symmetries
from ._helpers import CnScheme

_source = article(
    authors=["G.W. Tyler"],
    title="Numerical integration of functions of several variables",
    journal="Canad. J. Math.",
    volume="5",
    year="1953",
    pages="393-412",
    url="https://doi.org/10.4153/CJM-1953-044-1",
)


def tyler(n):
    d = {"0": [[frac(3 - n, 3)]], "a0": [[frac(1, 6)], [1]]}
    points, weights = expand_symmetries(d, n)
    return CnScheme("Tyler", n, weights, points, 3, _source)
