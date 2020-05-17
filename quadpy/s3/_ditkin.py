import sympy

from ..helpers import article, pm0, pm_roll, untangle
from ._helpers import S3Scheme

_source = article(
    authors=["V.A. Ditkin"],
    title="On certain approximate formulas for the calculation of triple integrals",
    journal="Doklady Akad. Nauk SSSR (N.S.)",
    number="62",
    year="1948",
    pages="445–447",
    note="Russian",
)

frac = sympy.Rational
sqrt = sympy.sqrt
pi = sympy.pi


def ditkin_1(alpha=0):
    B0 = frac(4, (alpha + 5) ** 2)
    B1 = frac((alpha + 3) * (alpha + 7), 12 * (alpha + 5) ** 2)

    r, s = [sqrt((alpha + 5) * (5 + i * sqrt(5)) / 10 / (alpha + 7)) for i in [+1, -1]]

    data = [
        (B0, [[0, 0, 0]]),
        (B1, pm_roll([r, s, 0])),
    ]

    points, weights = untangle(data)
    return S3Scheme("Ditkin 1", _source, 5, weights, points)


def ditkin_2():
    B0 = frac(4, 25)
    B1 = frac(21, 500)

    r, s = [sqrt((15 + i * 5 * sqrt(5)) / 42) for i in [+1, -1]]
    t = sqrt(frac(5, 21))

    data = [
        (B0, [[0, 0, 0]]),
        (B1, pm_roll([r, s, 0])),
        (B1, pm0([t, t, t])),
    ]

    points, weights = untangle(data)
    return S3Scheme("Ditkin 2", _source, 5, weights, points)


def ditkin_3():
    B0 = frac(16, 175)
    B1 = frac(81, 1400)
    B2 = frac(3, 280)

    sqrt5 = sqrt(5)
    r, s = [sqrt((5 + i * sqrt5) / 18) for i in [+1, -1]]
    t = sqrt(frac(1, 3))
    u, v = [sqrt((3 - i * sqrt5) / 6) for i in [+1, -1]]

    data = [
        (B0, [[0, 0, 0]]),
        (B1, pm_roll([r, s, 0])),
        (B2, pm_roll([u, v, 0])),
        (B2, pm0([t, t, t])),
    ]

    points, weights = untangle(data)
    return S3Scheme("Ditkin 3", _source, 7, weights, points)
