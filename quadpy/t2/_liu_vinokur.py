from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import T2Scheme, alpha, concat, gamma_delta, s3

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
    weights, points = s3(1)
    return T2Scheme("Liu-Vinokur 1", weights, points, 1, source)


def liu_vinokur_02():
    weights, points = alpha([frac(1, 3), 1])
    return T2Scheme("Liu-Vinokur 2", weights, points, 1, source)


def liu_vinokur_03():
    weights, points = alpha([frac(1, 3), -frac(1, 2)])
    return T2Scheme("Liu-Vinokur 3", weights, points, 2, source)


def liu_vinokur_04():
    weights, points = concat(s3(frac(3, 4)), alpha([frac(1, 12), 1]))
    return T2Scheme("Liu-Vinokur 4", weights, points, 2, source)


def liu_vinokur_05():
    weights, points = concat(
        s3(-frac(9, 16)),
        # ERR Incorrectly specified in the article as 25 (instead of 2/5).
        alpha([frac(25, 48), frac(2, 5)]),
    )
    return T2Scheme("Liu-Vinokur 5", weights, points, 3, source)


def liu_vinokur_06():
    sqrt21 = sqrt(21)
    weights, points = alpha(
        [(1 + sqrt21) / 120, 1], [(39 - sqrt21) / 120, (1 - sqrt21) / 10]
    )
    return T2Scheme("Liu-Vinokur 6", weights, points, 3, source)


def liu_vinokur_07():
    weights, points = concat(
        s3(frac(9, 20)), alpha([frac(1, 20), 1], [frac(2, 15), -frac(1, 2)])
    )
    return T2Scheme("Liu-Vinokur 7", weights, points, 3, source)


def liu_vinokur_08():
    sqrt10 = sqrt(10)
    sqrt_b = sqrt(950 - 220 * sqrt10)
    a1 = (-10 + 5 * sqrt10 + sqrt_b) / 30
    a2 = (-10 + 5 * sqrt10 - sqrt_b) / 30
    weights, points = alpha(
        [(5 * a2 - 2) / (60 * a1 ** 2 * (a2 - a1)), a1],
        [(5 * a1 - 2) / (60 * a2 ** 2 * (a1 - a2)), a2],
    )
    return T2Scheme("Liu-Vinokur 8", weights, points, 4, source)


def liu_vinokur_09():
    weights, points = concat(
        s3(frac(27, 80)),
        alpha([frac(8, 105), -frac(1, 2)], [frac(81, 560), frac(2, 3)]),
    )
    return T2Scheme("Liu-Vinokur 9", weights, points, 4, source)


def liu_vinokur_10():
    sqrt13 = sqrt(13)
    weights, points = alpha(
        [(11 - 1 * sqrt13) / 360, 1],
        [(80 - 16 * sqrt13) / 360, -frac(1, 2)],
        [(29 + 17 * sqrt13) / 360, (-1 + sqrt13) / 6],
    )
    return T2Scheme("Liu-Vinokur 10", weights, points, 4, source)


def liu_vinokur_11():
    sqrt3 = sqrt(3)
    weights, points = concat(
        s3(+frac(9, 20)),
        alpha([-frac(1, 60), 1]),
        gamma_delta([+frac(1, 10), (3 + sqrt3) / 6, (3 - sqrt3) / 6]),
    )
    return T2Scheme("Liu-Vinokur 11", weights, points, 4, source)


def liu_vinokur_12():
    sqrt15 = sqrt(15)
    weights, points = concat(
        s3(frac(9, 40)),
        alpha(
            [(155 - sqrt15) / 1200, (1 + sqrt15) / 7],
            [(155 + sqrt15) / 1200, (1 - sqrt15) / 7],
        ),
    )
    return T2Scheme("Liu-Vinokur 12", weights, points, 5, source)


def liu_vinokur_13():
    weights, points = concat(
        s3(frac(81, 320)),
        alpha(
            [frac(1, 90), 1],
            [frac(16, 225), -frac(1, 2)],
            [frac(2401, 14400), frac(4, 7)],
        ),
    )
    return T2Scheme("Liu-Vinokur 13", weights, points, 5, source)
