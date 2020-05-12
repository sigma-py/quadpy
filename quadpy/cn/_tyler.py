from sympy import Rational as frac

from ..helpers import article, fsd, untangle, z
from ._helpers import CnScheme

_citation = article(
    authors=["G.W. Tyler"],
    title="Numerical integration of functions of several variables",
    journal="Canad. J. Math.",
    volume="5",
    year="1953",
    pages="393-412",
    url="https://doi.org/10.4153/CJM-1953-044-1",
)


def tyler(n):
    data = [(frac(3 - n, 3), z(n)), (frac(1, 6), fsd(n, (1, 1)))]

    points, weights = untangle(data)
    reference_volume = 2 ** n
    weights *= reference_volume
    return CnScheme("Tyler", n, weights, points, 3, _citation)
