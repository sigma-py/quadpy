from sympy import Rational as frac

from ..helpers import article, expand_symmetries
from ._helpers import CnScheme

_source = article(
    authors=["G.M. Ewing"],
    title="On Approximate Cubature",
    journal="The American Mathematical Monthly",
    volume="48",
    number="2",
    month="feb",
    year="1941",
    pages="134-136",
    url="https://doi.org/10.2307/2303604",
)


def ewing(n):
    d = {"0": [[frac(2, 3)]], "a": [[frac(1, 3 * 2 ** n)], [1]]}
    points, weights = expand_symmetries(d, n)
    return CnScheme("Ewing", n, weights, points, 3, _source)
