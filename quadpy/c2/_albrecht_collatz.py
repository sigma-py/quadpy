from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, concat, pm, pm2, symm_r0, symm_s, zero

source = article(
    authors=["J. Albrecht", "L. Collatz"],
    title="Zur numerischen Auswertung mehrdimensionaler Integrale",
    journal="ZAMM",
    volume="38",
    number="1-2",
    year="1958",
    pages="1–15",
    url="https://doi.org/10.1002/zamm.19580380102",
)


def albrecht_collatz_1():
    weights, points = concat(
        zero(frac(5, 12)), symm_r0([frac(1, 8), 1]), symm_s([frac(1, 48), 1])
    )
    return C2Scheme("Albrecht-Collatz 1", weights, points, 3, source, 4.442e-16)


def albrecht_collatz_2():
    r = sqrt(frac(3, 5))
    s = sqrt(frac(1, 3))
    t = sqrt(frac(14, 15))
    weights, points = concat(
        zero(frac(2, 7)), pm([frac(5, 63), 0, t]), pm2([frac(5, 36), r, s])
    )
    return C2Scheme("Albrecht-Collatz 2", weights, points, 5, source, 4.627e-16)


def albrecht_collatz_3():
    r = sqrt(frac(7, 15))
    s, t = [sqrt((7 + i * sqrt(24)) / 15) for i in [+1, -1]]
    weights, points = concat(
        zero(frac(2, 7)),
        pm([frac(25, 168), r, r], [frac(5, 48), +s, -t], [frac(5, 48), +t, -s]),
    )
    return C2Scheme("Albrecht-Collatz 3", weights, points, 5, source, 4.442e-16)


def albrecht_collatz_4():
    weights, points = concat(
        zero(frac(2, 45)),
        symm_r0([frac(2, 45), 1]),
        symm_s([frac(1, 60), 1], [frac(8, 45), frac(1, 2)]),
    )
    return C2Scheme("Albrecht-Collatz 4", weights, points, 5, source, 8.883e-16)
