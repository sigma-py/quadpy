import numpy as np
import sympy
from mpmath import mp

from ..helpers import article
from ._helpers import S2Scheme, register

_source = article(
    authors=["J. Albrecht"],
    title="Formeln zur numerischen Integration über Kreisbereiche",
    journal="ZAMM",
    volume="40",
    number="10-11",
    year="1960",
    pages="514–517",
    url="https://doi.org/10.1002/zamm.19600401014",
)

frac = sympy.Rational
pi = sympy.pi
cos = np.vectorize(sympy.cos)
sin = np.vectorize(sympy.sin)
sqrt = np.vectorize(sympy.sqrt)
pm_ = np.array([+1, -1])
roots = mp.polyroots
linear_solve = mp.lu_solve


def albrecht_1():
    # Equals Albrecht-Collatz, Lether(2)
    d = {"d4_aa": [[frac(1, 4)], [frac(1, 2)]]}
    return S2Scheme("Albrecht 1", d, 3, _source)


def albrecht_2():
    d = {"zero2": [[frac(1, 4)]], "d6.1": [[frac(1, 8)], [sqrt(frac(2, 3))]]}
    return S2Scheme("Albrecht 2", d, 5, _source)


def albrecht_3():
    sqrt29 = sqrt(29)
    a1, a2 = (551 + pm_ * 41 * sqrt29) / 6264
    rho1, rho2 = sqrt((27 - pm_ * 3 * sqrt29) / 52)

    d = {"d4.1": [[frac(2, 27)], [sqrt(frac(3, 4))]], "d4.0": [[a1, a2], [rho1, rho2]]}
    return S2Scheme("Albrecht 3", d, 7, _source)


def albrecht_4():
    sqrt111 = sqrt(111)
    rho1, rho2 = sqrt((96 - pm_ * 4 * sqrt(111)) / 155)

    B0 = frac(251, 2304)
    B1, B2 = (110297 + pm_ * 5713 * sqrt111) / 2045952
    C = frac(125, 3072)

    d = {
        "zero2": [[B0]],
        "d6.0": [[B1, B2], [rho1, rho2]],
        "d6.1": [[C], [sqrt(frac(4, 5))]],
    }
    return S2Scheme("Albrecht 4", d, 9, _source)


def albrecht_5():
    # The values are solutions of
    # 6317094x^3 - 10022245*x^2 + 4149900*x - 336375 = 0
    sigma2 = roots([6317094, -10022245, 4149900, -336375])
    A = np.vander(sigma2, increasing=True).T
    b = np.array([frac(168899, 1350000), frac(7661, 180000), frac(71, 3000)])
    B = linear_solve(A, b)

    sqrt19 = sqrt(19)

    # ERR Stroud incorrectly lists sqrt(10) for s1.
    s1, s2 = sqrt((125 - pm_ * 10 * sqrt19) / 366)

    # ERR Stroud incorrectly lists 749489_3_.0 instead of 749489_2_.0
    C1, C2 = (7494892 + pm_ * 1053263 * sqrt19) / 205200000
    D = frac(81, 3125)

    u = sqrt(frac(5, 6)) * cos(pi / 8)
    v = sqrt(frac(5, 6)) * sin(pi / 8)

    d = {
        "d4_a0": [B, [sqrt(sigma2[0]), sqrt(sigma2[1]), sqrt(sigma2[2])]],
        "d4_aa": [[C1, C2], [s1, s2]],
        "d4_ab": [[D], [u], [v]],
    }
    return S2Scheme("Albrecht 5", d, 11, _source)


def albrecht_6():
    # The values are solutions of
    # 11025*x^3 - 19020*x^2 + 9370*x - 1212 = 0
    sigma2 = roots([11025, -19020, 9370, -1212])
    A = np.vander(sigma2, increasing=True).T
    b = np.array([frac(1432433, 18849024), frac(1075, 31104), frac(521, 25920)])
    B = linear_solve(A, b)

    B0 = frac(2615, 43632)
    C = frac(16807, 933120)

    d = {
        "zero2": [[B0]],
        "d10.0": [B, sqrt(sigma2)],
        "d10.1": [[C], [sqrt(frac(6, 7))]],
    }
    return S2Scheme("Albrecht 6", d, 13, _source)


def albrecht_7():
    sqrt21 = sqrt(21)
    wt1, wt2 = (4998 + pm_ * 343 * sqrt21) / 253125
    tau1, tau2 = sqrt((21 - pm_ * sqrt21) / 28)

    # The values are solutions of
    # 4960228*x^4 - 10267740*x^3 + 6746490*x^2 - 1476540*x + 70425 = 0
    sigma2 = roots([4960228, -10267740, 6746490, -1476540, 70425])
    A = np.vander(sigma2, increasing=True).T
    b = np.array(
        [frac(57719, 675000), frac(9427, 270000), frac(193, 9000), frac(113, 7200)]
    )
    ws = linear_solve(A, b)

    d = {
        "d8.0": [ws, sqrt(sigma2)],
        "d8.1": [[wt1, wt2], [tau1, tau2]],
    }
    return S2Scheme("Albrecht 7", d, 15, _source)


def albrecht_8():
    m0 = frac(496439663, 13349499975)

    sqrt7 = sqrt(7)
    wt1, wt2 = (125504 + pm_ * 16054 * sqrt7) / 8751645
    tau1, tau2 = sqrt((14 - pm_ * sqrt7) / 18)

    # The values are solutions of
    # 160901628*x^4 - 364759920*x^3 + 274856190*x^2 - 76570340*x
    # + 6054195 = 0
    sigma2 = roots([160901628, -364759920, 274856190, -76570340, 6054195])
    A = np.vander(sigma2, increasing=True).T
    b = np.array(
        [
            frac(121827491812, 1802182496625),
            frac(48541, 1666980),
            frac(977, 55566),
            frac(671, 52920),
        ]
    )
    ws = linear_solve(A, b)

    d = {
        "zero2": [[m0]],
        "d10.0": [ws, sqrt(sigma2)],
        "d10.1": [[wt1, wt2], [tau1, tau2]],
    }
    return S2Scheme("Albrecht 8", d, 17, _source)


register(
    [
        albrecht_1,
        albrecht_2,
        albrecht_3,
        albrecht_4,
        albrecht_5,
        albrecht_6,
        albrecht_7,
        albrecht_8,
    ]
)
