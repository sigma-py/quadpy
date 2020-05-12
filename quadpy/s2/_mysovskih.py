import numpy
import sympy

from ..helpers import article, fs_array, fsd, untangle, z
from ._helpers import S2Scheme

_citation = article(
    authors=["I.P. Mysovskikh"],
    title="On the construction of cubature formulas for the simplest regions",
    journal="Z. Vychisl. Mat. i. Mat. Fiz.",
    number="4",
    pages="3-14",
    year="1964",
)

frac = sympy.Rational
pi = sympy.pi
sqrt = numpy.vectorize(sympy.sqrt)
cos = numpy.vectorize(sympy.cos)
sin = numpy.vectorize(sympy.sin)
pm_ = numpy.array([+1, -1])


def mysovskih_1(alpha=0):
    b = sqrt(frac(alpha + 4, alpha + 6))

    a = 2 * numpy.arange(5) * pi / 5
    x = b * numpy.array([cos(a), sin(a)]).T

    B0 = frac(4, (alpha + 4) ** 2)
    B1 = frac((alpha + 2) * (alpha + 6), 5 * (alpha + 4) ** 2)

    data = [(B0, z(2)), (B1, x)]

    points, weights = untangle(data)
    weights *= pi
    return S2Scheme("Mysovskih 1", weights, points, 4, _citation)


def mysovskih_2():
    sqrt10 = sqrt(10)
    sqrt601 = sqrt(601)

    B1, B3 = (857 * sqrt601 + pm_ * 12707) / 20736 / sqrt601
    B2 = frac(125, 3456)
    B4, B5 = (340 + pm_ * 25 * sqrt10) / 10368

    r1, r3 = sqrt((31 - pm_ * sqrt601) / 60)
    r2 = sqrt(frac(3, 5))
    r4, r5 = sqrt((10 - pm_ * sqrt10) / 20)

    s4, s5 = sqrt((10 - pm_ * sqrt10) / 60)

    data = [
        (B1, fsd(2, (r1, 1))),
        (B2, fsd(2, (r2, 1))),
        (B3, fsd(2, (r3, 1))),
        (B4, fs_array([r4, s4])),
        (B5, fs_array([r5, s5])),
    ]

    points, weights = untangle(data)
    weights *= pi
    return S2Scheme("Mysovskih 2", weights, points, 11, _citation)


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

    a = (2 * numpy.arange(8) + 1) * pi / 8
    xa = numpy.array([cos(a), sin(a)]).T

    a = (2 * numpy.arange(4) + 1) * pi / 4
    xb = numpy.array([cos(a), sin(a)]).T

    a = numpy.arange(1, 5) * pi / 2
    xc = numpy.array([cos(a), sin(a)]).T

    data = [
        (A1, rho1 * xa),
        (A2, rho2 * xa),
        (B1, sigma1 * xb),
        (B2, sigma2 * xb),
        (B3, sigma3 * xb),
        (C1, tau1 * xc),
        (C2, tau2 * xc),
        (C3, tau3 * xc),
        (C4, tau4 * xc),
    ]

    points, weights = untangle(data)
    weights *= pi
    return S2Scheme("Mysovskih 3", weights, points, 15, _citation)
