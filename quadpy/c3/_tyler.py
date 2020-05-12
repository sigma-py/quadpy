from sympy import Rational as frac

from ..helpers import article, untangle
from ._helpers import HexahedronScheme, fs_r00, pm_rrr, z

citation = article(
    authors=["G.W. Tyler"],
    title="Numerical integration of functions of several variables",
    journal="Canad. J. Math.",
    volume="5",
    year="1953",
    pages="393-412",
    url="https://doi.org/10.4153/CJM-1953-044-1",
)


def tyler_1():
    data = [(frac(1, 6), fs_r00(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Tyler 1", weights, points, 3, citation)


def tyler_2():
    data = [
        (-frac(62, 45), z()),
        (frac(16, 45), fs_r00(frac(1, 2))),
        (frac(1, 45), fs_r00(1)),
        (frac(1, 72), pm_rrr(1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Tyler 2", weights, points, 5, citation)
