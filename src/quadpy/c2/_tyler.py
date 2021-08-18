from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, register

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
    d = {
        "zero2": [[-frac(28, 45)]],
        "d4_aa": [[frac(1, 36)], [1]],
        "d4_a0": [[frac(1, 45), frac(16, 45)], [1, frac(1, 2)]],
    }
    return C2Scheme("Tyler 1", d, 5, source)


def tyler_2():
    r = sqrt(frac(6, 7))
    s, t = (sqrt((114 - i * 3 * sqrt(583)) / 287) for i in [+1, -1])
    B1 = frac(49, 810)
    B2, B3 = ((178981 + i * 2769 * sqrt(583)) / 1888920 for i in [+1, -1])
    d = {"d4_a0": [[B1], [r]], "d4_aa": [[B2, B3], [s, t]]}
    return C2Scheme("Tyler 2", d, 7, source, 3.331e-16)


def tyler_3():
    d = {
        "zero2": [[frac(449, 315)]],
        "d4_a0": [
            [frac(37, 1260), frac(3, 28), -frac(69, 140)],
            [1, frac(2, 3), frac(1, 3)],
        ],
        "d4_aa": [[frac(7, 540), frac(32, 135)], [1, frac(1, 2)]],
    }
    return C2Scheme("Tyler 3", d, 7, source)


register([tyler_1, tyler_2, tyler_3])
