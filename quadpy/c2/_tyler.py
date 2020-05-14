from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, concat, symm_r0, symm_s, zero

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
    weights, points = concat(
        zero(-frac(28, 45)),
        symm_s([frac(1, 36), 1]),
        symm_r0([frac(1, 45), 1], [frac(16, 45), frac(1, 2)]),
    )
    return C2Scheme("Tyler 1", weights, points, 5, source)


def tyler_2():
    r = sqrt(frac(6, 7))
    s, t = [sqrt((114 - i * 3 * sqrt(583)) / 287) for i in [+1, -1]]
    B1 = frac(49, 810)
    B2, B3 = [(178981 + i * 2769 * sqrt(583)) / 1888920 for i in [+1, -1]]
    weights, points = concat(symm_r0([B1, r]), symm_s([B2, s], [B3, t]))
    return C2Scheme("Tyler 2", weights, points, 7, source)


def tyler_3():
    weights, points = concat(
        zero(frac(449, 315)),
        symm_r0(
            [frac(37, 1260), 1], [frac(3, 28), frac(2, 3)], [-frac(69, 140), frac(1, 3)]
        ),
        symm_s([frac(7, 540), 1], [frac(32, 135), frac(1, 2)]),
    )
    return C2Scheme("Tyler 3", weights, points, 7, source)
