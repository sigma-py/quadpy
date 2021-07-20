import pathlib

from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import T3Scheme, _read, register

source = article(
    authors=["P. Keast"],
    title="Moderate degree tetrahedral quadrature formulas",
    journal="Computer Methods in Applied Mechanics and Engineering",
    volume="55",
    number="3",
    pages="339-348",
    month="may",
    year="1986",
    url="https://doi.org/10.1016/0045-7825%2886%2990059-9",
)

# https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
this_dir = pathlib.Path(__file__).resolve().parent


def keast_0():
    # Does not appear in Keast's article. TODO remove
    degree = 1
    d = {"s4": [[1]]}
    return T3Scheme("Keast 0", d, degree, source)


def keast_1():
    # Does not appear in Keast's article.
    degree = 2
    d = {"s31": [[frac(1, 4)], [frac(1, 4) - sqrt(5) / 20]]}
    return T3Scheme("Keast 1", d, degree, source)


def keast_2():
    # Does not appear in Keast's article.
    degree = 3
    d = {"s4": [[-frac(4, 5)]], "s31": [[frac(9, 20)], [frac(1, 6)]]}
    return T3Scheme("Keast 2", d, degree, source)


def keast_3():
    # Does not appear in Keast's article.
    degree = 3
    d = {
        "s31": [[0.2177650698804054], [0.1438564719343852]],
        "s22": [[0.0214899534130631], [0.5]],
    }
    return T3Scheme("Keast 3", d, degree, source)


def keast_4():
    degree = 4
    d = {
        "s4": [[-frac(148, 1875)]],
        "s31": [[frac(343, 7500)], [frac(1, 14)]],
        "s22": [[frac(56, 375)], [frac(1, 4) + sqrt(frac(5, 14)) / 4]],
    }
    return T3Scheme("Keast 4", d, degree, source)


def keast_5():
    degree = 4
    d = {
        "s22": [[2 / 105], [0.5]],
        "s31": [
            [0.0885898247429807, 0.1328387466855907],
            [0.1005267652252045, 0.3143728734931922],
        ],
    }
    return T3Scheme("Keast 5", d, degree, source)


def keast_6():
    degree = 5
    d = {
        "s4": [[frac(6544, 36015)]],
        "s31": [[frac(81, 2240), frac(161051, 2304960)], [frac(1, 3), frac(1, 11)]],
        "s22": [[frac(338, 5145)], [frac(1, 4) - sqrt(91) / 52]],
    }
    return T3Scheme("Keast 6", d, degree, source)


def keast_7():
    return _read(this_dir / "keast_7.json", source)


def keast_8():
    return _read(this_dir / "keast_8.json", source)


def keast_9():
    return _read(this_dir / "keast_9.json", source)


register(
    [
        keast_0,
        keast_1,
        keast_2,
        keast_3,
        keast_4,
        keast_5,
        keast_6,
        keast_7,
        keast_8,
        keast_9,
    ]
)
