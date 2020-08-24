from sympy import Rational as frac

from ..helpers import article
from ._helpers import C3Scheme, register

source = article(
    authors=["G.W. Tyler"],
    title="Numerical integration of functions of several variables",
    journal="Canad. J. Math.",
    volume="5",
    year="1953",
    pages="393-412",
    url="https://doi.org/10.4153/CJM-1953-044-1",
)


def tyler_1():
    d = {"symm_r00": [[frac(1, 6)], [1]]}
    return C3Scheme("Tyler 1", d, 3, source)


def tyler_2():
    d = {
        "zero3": [[-frac(62, 45)]],
        "symm_r00": [[frac(16, 45), frac(1, 45)], [frac(1, 2), 1]],
        "symm_rrr": [[frac(1, 72)], [1]],
    }
    return C3Scheme("Tyler 2", d, 5, source)


register([tyler_1, tyler_2])
