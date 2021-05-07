import numpy as np
import sympy

from ..helpers import article
from ._helpers import S2Scheme, register

_source = article(
    authors=["I.P. Mysovskikh"],
    title="On the construction of cubature formulas for the simplest regions",
    journal="Z. Vychisl. Mat. i. Mat. Fiz.",
    number="4",
    pages="3-14",
    year="1964",
)

frac = sympy.Rational
pi = sympy.pi
sqrt = np.vectorize(sympy.sqrt)
cos = np.vectorize(sympy.cos)
sin = np.vectorize(sympy.sin)
pm_ = np.array([+1, -1])


def mysovskih_1(alpha=0):
    b = sqrt(frac(alpha + 4, alpha + 6))

    B0 = frac(4, (alpha + 4) ** 2)
    B1 = frac((alpha + 2) * (alpha + 6), 5 * (alpha + 4) ** 2)

    d = {
        "zero2": [[B0]],
        "d5.0": [[B1], [b]],
    }
    return S2Scheme("Mysovskih 1", d, 4, _source)


def mysovskih_2():
    sqrt10 = sqrt(10)
    sqrt601 = sqrt(601)

    B1, B3 = (857 * sqrt601 + pm_ * 12707) / 20736 / sqrt601
    B2 = frac(125, 3456)
    B4, B5 = (340 + pm_ * 25 * sqrt10) / 10368

    r1, r3 = sqrt((31 - pm_ * sqrt601) / 60)
    r2 = sqrt(frac(3, 5))
    r4, r5 = sqrt((10 - pm_ * sqrt10) / 20)

    c4, s5 = sqrt((10 - pm_ * sqrt10) / 60)

    d = {"d4_a0": [[B1, B2, B3], [r1, r2, r3]], "d4_ab": [[B4, B5], [r4, r5], [c4, s5]]}
    return S2Scheme("Mysovskih 2", d, 11, _source)


def mysovskih_3():
    sqrt21 = sqrt(21)
    sqrt1401 = sqrt(1401)

    # ERR Stroud lists 5096 instead of 4998 here
    A1, A2 = (4998 + pm_ * 343 * sqrt21) / 253125

    # ERR Stroud is missing the +- here
    B1, B3 = (1055603 * sqrt1401 + pm_ * 26076047) / 43200000 / sqrt1401

    B2 = frac(16807, 800000)

    rho1, rho2 = sqrt((21 - pm_ * sqrt21) / 28)
    sigma1, sigma3 = sqrt((69 - pm_ * sqrt1401) / 112)
    sigma2 = sqrt(frac(5, 7))

    tau1 = 0.252863797091230
    tau2 = 0.577728928444823
    tau3 = 0.873836956644882
    tau4 = 0.989746802511491

    C1 = 0.398811120280412e-1
    C2 = 0.348550570365141e-1
    C3 = 0.210840370156484e-1
    C4 = 0.531979391979623e-2

    d = {
        "d8.1": [[A1, A2], [rho1, rho2]],
        "d4.1": [[B1, B2, B3], [sigma1, sigma2, sigma3]],
        "d4.0": [[C1, C2, C3, C4], [tau1, tau2, tau3, tau4]],
    }
    return S2Scheme("Mysovskih 3", d, 15, _source)


register([mysovskih_1, mysovskih_2, mysovskih_3])
