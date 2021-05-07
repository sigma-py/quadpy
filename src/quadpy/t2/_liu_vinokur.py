from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._classical import centroid as liu_vinokur_01
from ._classical import seven_point as liu_vinokur_07
from ._classical import vertex as liu_vinokur_02
from ._helpers import T2Scheme, register

source = article(
    authors=["Y. Liu", "M. Vinokur"],
    title="Exact Integrations of Polynomials and Symmetric Quadrature Formulas over Arbitrary Polyhedral Grids",
    journal="Journal of Computational Physics",
    volume="140",
    pages="122–147",
    year="1998",
    url="https://doi.org/10.1006/jcph.1998.5884",
)


def liu_vinokur_03():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 2)]]}
    return T2Scheme("Liu-Vinokur 3", d, 2, source)


def liu_vinokur_04():
    d = {"centroid": [[frac(3, 4)]], "vertex": [[frac(1, 12)]]}
    return T2Scheme("Liu-Vinokur 4", d, 2, source)


def liu_vinokur_05():
    # ERR Incorrectly specified in the article as 25 (instead of 2/5).
    # alpha = frac(2, 5)
    # b = (1 - alpha) / 3
    d = {"centroid": [[-frac(9, 16)]], "d3_aa": [[frac(25, 48)], [frac(1, 5)]]}
    return T2Scheme("Liu-Vinokur 5", d, 3, source)


def liu_vinokur_06():
    sqrt21 = sqrt(21)
    alpha1 = (1 - sqrt21) / 10
    b1 = (1 - alpha1) / 3
    d = {"vertex": [[(1 + sqrt21) / 120]], "d3_aa": [[(39 - sqrt21) / 120], [b1]]}
    return T2Scheme("Liu-Vinokur 6", d, 3, source)


def liu_vinokur_08():
    sqrt10 = sqrt(10)
    sqrt_b = sqrt(950 - 220 * sqrt10)
    a1 = (-10 + 5 * sqrt10 + sqrt_b) / 30
    a2 = (-10 + 5 * sqrt10 - sqrt_b) / 30
    b1 = (1 - a1) / 3
    b2 = (1 - a2) / 3
    d = {
        "d3_aa": [
            [
                (5 * a2 - 2) / (60 * a1 ** 2 * (a2 - a1)),
                (5 * a1 - 2) / (60 * a2 ** 2 * (a1 - a2)),
            ],
            [b1, b2],
        ]
    }
    return T2Scheme("Liu-Vinokur 8", d, 4, source)


def liu_vinokur_09():
    alpha0 = -frac(1, 2)
    alpha1 = frac(2, 3)
    b0 = (1 - alpha0) / 3
    b1 = (1 - alpha1) / 3
    d = {
        "centroid": [[frac(27, 80)]],
        "d3_aa": [[frac(8, 105), frac(81, 560)], [b0, b1]],
    }
    return T2Scheme("Liu-Vinokur 9", d, 4, source)


def liu_vinokur_10():
    sqrt13 = sqrt(13)

    alpha0 = 1
    alpha1 = -frac(1, 2)
    alpha2 = (-1 + sqrt13) / 6
    b0 = (1 - alpha0) / 3
    b1 = (1 - alpha1) / 3
    b2 = (1 - alpha2) / 3
    d = {
        "d3_aa": [
            [
                (11 - 1 * sqrt13) / 360,
                (80 - 16 * sqrt13) / 360,
                (29 + 17 * sqrt13) / 360,
            ],
            [b0, b1, b2],
        ]
    }
    return T2Scheme("Liu-Vinokur 10", d, 4, source)


def liu_vinokur_11():
    sqrt3 = sqrt(3)

    gamma = (3 + sqrt3) / 6
    delta = (3 - sqrt3) / 6
    a = (1 + 2 * gamma - delta) / 3
    b = (1 + 2 * delta - gamma) / 3
    # c = (1 - gamma - delta) / 3

    d = {
        "centroid": [[frac(9, 20)]],
        "vertex": [[-frac(1, 60)]],
        "d3_ab": [[frac(1, 10)], [a], [b]],
    }
    return T2Scheme("Liu-Vinokur 11", d, 4, source)


def liu_vinokur_12():
    sqrt15 = sqrt(15)

    a0 = (1 + sqrt15) / 7
    a1 = (1 - sqrt15) / 7
    b0 = (1 - a0) / 3
    b1 = (1 - a1) / 3
    d = {
        "centroid": [[frac(9, 40)]],
        "d3_aa": [[(155 - sqrt15) / 1200, (155 + sqrt15) / 1200], [b0, b1]],
    }
    return T2Scheme("Liu-Vinokur 12", d, 5, source)


def liu_vinokur_13():
    d = {
        "centroid": [[frac(81, 320)]],
        "vertex": [[frac(1, 90)]],
        "d3_aa": [[frac(16, 225), frac(2401, 14400)], [frac(1, 2), frac(1, 7)]],
    }
    return T2Scheme("Liu-Vinokur 13", d, 5, source)


register(
    [
        liu_vinokur_01,
        liu_vinokur_02,
        liu_vinokur_03,
        liu_vinokur_04,
        liu_vinokur_05,
        liu_vinokur_06,
        liu_vinokur_07,
        liu_vinokur_08,
        liu_vinokur_09,
        liu_vinokur_10,
        liu_vinokur_11,
        liu_vinokur_12,
        liu_vinokur_13,
    ]
)
