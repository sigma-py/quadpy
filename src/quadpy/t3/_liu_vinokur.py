from sympy import Rational as frac
from sympy import acos, cos, pi, sqrt

from ..helpers import article
from ._helpers import T3Scheme, register

source = article(
    authors=["Y. Liu", "M. Vinokur"],
    title="Exact Integrations of Polynomials and Symmetric Quadrature Formulas over Arbitrary Polyhedral Grids",
    journal="Journal of Computational Physics",
    volume="140",
    pages="122–147",
    year="1998",
    url="https://doi.org/10.1006/jcph.1998.5884",
)


def liu_vinokur_01():
    d = {"s4": [[1]]}
    degree = 1
    return T3Scheme("Liu-Vinokur 1", d, degree, source)


def liu_vinokur_02():
    d = {"s31": [[frac(1, 4)], [0]]}
    degree = 1
    return T3Scheme("Liu-Vinokur 2", d, degree, source)


def liu_vinokur_03():
    d = {"s31": [[frac(1, 4)], [frac(1, 4) - sqrt(5) / 20]]}
    degree = 2
    return T3Scheme("Liu-Vinokur 3", d, degree, source)


def liu_vinokur_04():
    d = {"s4": [[frac(4, 5)]], "s31": [[frac(1, 20)], [0]]}
    degree = 2
    return T3Scheme("Liu-Vinokur 4", d, degree, source)


def liu_vinokur_05():
    d = {"s4": [[-frac(4, 5)]], "s31": [[frac(9, 20)], [frac(1, 6)]]}
    degree = 3
    return T3Scheme("Liu-Vinokur 5", d, degree, source)


def liu_vinokur_06():
    d = {"s31": [[frac(1, 40), frac(9, 40)], [0, frac(1, 3)]]}
    degree = 3
    return T3Scheme("Liu-Vinokur 6", d, degree, source)


def liu_vinokur_07():
    b = (1 + 2 * sqrt(70) / 28) / 4
    d = {
        "s4": [[-frac(148, 1875)]],
        "s31": [[frac(343, 7500)], [frac(1, 14)]],
        "s22": [[frac(56, 375)], [b]],
    }
    degree = 4
    return T3Scheme("Liu-Vinokur 7", d, degree, source)


def liu_vinokur_08():
    alpha1 = (+sqrt(65944 - 19446 * sqrt(11)) + 51 * sqrt(11) - 154) / 89
    alpha2 = (-sqrt(65944 - 19446 * sqrt(11)) + 51 * sqrt(11) - 154) / 89
    d = {
        "s31": [
            [
                (17 * alpha2 - 7) / (420 * alpha1 ** 2 * (alpha2 - alpha1)),
                (17 * alpha1 - 7) / (420 * alpha2 ** 2 * (alpha1 - alpha2)),
            ],
            [(1 - alpha1) / 4, (1 - alpha2) / 4],
        ],
        "s22": [[frac(2, 105)], [0]],
    }
    degree = 4
    return T3Scheme("Liu-Vinokur 8", d, degree, source)


def liu_vinokur_09():
    d = {
        "s4": [[-frac(32, 15)]],
        "s31": [[frac(3, 280), frac(125, 168)], [0, frac(1, 5)]],
        "s22": [[frac(2, 105)], [0]],
    }
    degree = 4
    return T3Scheme("Liu-Vinokur 9", d, degree, source)


def liu_vinokur_10():
    gamma1 = (2 + sqrt(2)) / 4
    gamma2 = (2 - sqrt(2)) / 4
    a = (1 - gamma1 - gamma2) / 4
    b = (1 + 3 * gamma1 - gamma2) / 4
    # c = (1 + 3 * delta - gamma) / 4
    d = {
        "s4": [[frac(32, 105)]],
        "s31": [[-frac(31, 840), frac(27, 280)], [0, frac(1, 3)]],
        "s211": [[frac(4, 105)], [a], [b]],
    }
    degree = 4
    return T3Scheme("Liu-Vinokur 10", d, degree, source)


def liu_vinokur_11():
    d = {
        "s31": [
            [
                (11 - 4 * sqrt(2)) / 840,
                (243 - 108 * sqrt(2)) / 1960,
                (62 + 44 * sqrt(2)) / 735,
            ],
            [0, frac(1, 3), (2 - sqrt(2)) / 4],
        ],
        "s22": [[frac(2, 105)], [0]],
    }
    degree = 4
    return T3Scheme("Liu-Vinokur 11", d, degree, source)


def liu_vinokur_12():
    lmbda = frac(4, 27) * (
        4 * sqrt(79) * cos((acos(67 * sqrt(79) / 24964) + 2 * pi) / 3) + 71
    )
    alpha1 = (+sqrt(9 * lmbda ** 2 - 248 * lmbda + 1680) + 28 - 3 * lmbda) / (
        112 - 10 * lmbda
    )
    alpha2 = (-sqrt(9 * lmbda ** 2 - 248 * lmbda + 1680) + 28 - 3 * lmbda) / (
        112 - 10 * lmbda
    )
    w1 = ((21 - lmbda) * alpha2 - 7) / (420 * alpha1 ** 2 * (alpha2 - alpha1))
    w2 = ((21 - lmbda) * alpha1 - 7) / (420 * alpha2 ** 2 * (alpha1 - alpha2))
    d = {
        "s31": [[w1, w2], [(1 - alpha1) / 4, (1 - alpha2) / 4]],
        "s22": [[lmbda ** 2 / 840], [(1 - 2 / sqrt(lmbda)) / 4]],
    }
    degree = 5
    return T3Scheme("Liu-Vinokur 12", d, degree, source)


def liu_vinokur_13():
    d = {
        "s4": [[-frac(16, 21)]],
        "s31": [
            [(2249 - 391 * sqrt(13)) / 10920, (2249 + 391 * sqrt(13)) / 10920],
            [(1 - (2 + sqrt(13)) / 9) / 4, (1 - (2 - sqrt(13)) / 9) / 4],
        ],
        "s22": [[frac(2, 105)], [0]],
    }
    degree = 5
    return T3Scheme("Liu-Vinokur 13", d, degree, source)


def liu_vinokur_14():
    d = {
        "s4": [[frac(16, 105)]],
        "s31": [
            [frac(1, 280), frac(81, 1400), frac(64, 525)],
            [0, frac(1, 3), frac(1, 8)],
        ],
        "s22": [[frac(2, 105)], [0]],
    }
    degree = 5
    return T3Scheme("Liu-Vinokur 14", d, degree, source)


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
        liu_vinokur_14,
    ]
)
