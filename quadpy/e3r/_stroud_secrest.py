from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, pm, pm_roll, untangle
from ._helpers import E3rScheme

source = article(
    authors=["A.H. Stroud", "D. Secrest"],
    title="Approximate integration formulas for certain spherically symmetric regions",
    journal="Math. Comp.",
    volume="17",
    year="1963",
    pages="105-135",
    url="https://doi.org/10.1090/S0025-5718-1963-0161473-0",
)


def stroud_secrest_07():
    nu, xi = [sqrt(15 - p_m * 3 * sqrt(5)) for p_m in [+1, -1]]
    A = frac(3, 5)
    B = frac(1, 30)

    data = [(A, [[0, 0, 0]]), (B, pm_roll([xi, nu, 0]))]

    points, weights = untangle(data)
    return E3rScheme("Stroud-Secrest VII", weights, points, 5, source)


def stroud_secrest_08():
    nu = sqrt(30)
    eta = sqrt(10)
    A = frac(3, 5)
    B = frac(2, 75)
    C = frac(3, 100)

    data = [(A, [[0, 0, 0]]), (B, pm_roll([nu, 0, 0])), (C, pm([eta, eta, eta]))]
    points, weights = untangle(data)
    return E3rScheme("Stroud-Secrest VIII", weights, points, 5, source)


def stroud_secrest_09():
    eta = sqrt(10)
    xi, nu = [sqrt(15 - p_m * 5 * sqrt(5)) for p_m in [+1, -1]]
    A = frac(3, 5)
    B = frac(1, 50)

    data = [(A, [[0, 0, 0]]), (B, pm([eta, eta, eta])), (B, pm_roll([xi, nu, 0]))]
    points, weights = untangle(data)
    return E3rScheme("Stroud-Secrest IX", weights, points, 5, source)


def stroud_secrest_10():
    sqrt130 = sqrt(130)

    nu = sqrt((720 - 24 * sqrt130) / 11)
    xi = sqrt(288 + 24 * sqrt130)
    eta = sqrt((-216 + 24 * sqrt130) / 7)
    A = (5175 - 13 * sqrt130) / 8820
    B = (3870 + 283 * sqrt130) / 493920
    C = (3204 - 281 * sqrt130) / 197568
    # ERR in Stroud's book: 917568 vs. 197568
    D = (4239 + 373 * sqrt130) / 197568

    data = [
        (A, [[0, 0, 0]]),
        (B, fsd(3, (nu, 1))),
        (C, fsd(3, (xi, 2))),
        (D, pm([eta, eta, eta])),
    ]

    points, weights = untangle(data)
    return E3rScheme("Stroud-Secrest X", weights, points, 7, source)


def stroud_secrest_11():
    sqrt5 = sqrt(5)
    sqrt39 = sqrt(39)
    sqrt195 = sqrt(195)

    nu, xi = [
        sqrt(-50 + p_m * 10 * sqrt5 + 10 * sqrt39 - p_m * 2 * sqrt195)
        for p_m in [+1, -1]
    ]
    eta = sqrt(36 + 4 * sqrt39)
    mu, lmbda = [
        sqrt(54 + p_m * 18 * sqrt5 + 6 * sqrt39 + p_m * 2 * sqrt195) for p_m in [+1, -1]
    ]
    A = (1725 - 26 * sqrt39) / 2940
    B = (1065 + 171 * sqrt39) / 54880
    C = (297 - 47 * sqrt39) / 32928

    data = [
        (A, [[0, 0, 0]]),
        (B, pm_roll([xi, nu, 0])),
        (C, pm([eta, eta, eta])),
        (C, pm_roll([lmbda, mu, 0])),
    ]

    points, weights = untangle(data)
    return E3rScheme("Stroud-Secrest XI", weights, points, 7, source, 1.487e-11)
