from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import U3Scheme, register

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
    r, s = (sqrt((5 + t * sqrt(5)) / 10) for t in [+1, -1])
    data = {"rs0": [frac(1, 12), r, s]}
    return U3Scheme("Albrecht-Collatz 1", data, 5, source)


def albrecht_collatz_2():
    data = {"a1": [frac(8, 120)], "a3": [frac(9, 120)]}
    return U3Scheme("Albrecht-Collatz 2", data, 5, source)


def albrecht_collatz_3():
    data = {"a1": [frac(1, 30)], "a2": [frac(2, 30)]}
    return U3Scheme("Albrecht-Collatz 3", data, 5, source)


def albrecht_collatz_4():
    r, s = (sqrt((3 + t * sqrt(5)) / 6) for t in [+1, -1])
    data = {"rs0": [frac(1, 20), r, s], "a3": [frac(1, 20)]}
    return U3Scheme("Albrecht-Collatz 4", data, 5, source)


def albrecht_collatz_5():
    data = {
        "a1": [frac(40, 840)],
        "a2": [frac(32, 840)],
        "a3": [frac(27, 840)],
    }
    return U3Scheme("Albrecht-Collatz 5", data, 7, source)


register(
    [
        albrecht_collatz_1,
        albrecht_collatz_2,
        albrecht_collatz_3,
        albrecht_collatz_4,
        albrecht_collatz_5,
    ]
)
