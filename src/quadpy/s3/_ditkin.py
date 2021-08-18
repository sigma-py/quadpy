import sympy

from ..helpers import article
from ._helpers import S3Scheme, register

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

    r, s = (sqrt((alpha + 5) * (5 + i * sqrt(5)) / 10 / (alpha + 7)) for i in [+1, -1])

    d = {"zero3": [[B0]], "symm_rs0_roll": [[B1], [r], [s]]}
    return S3Scheme("Ditkin 1", d, 5, _source)


def ditkin_2():
    B0 = frac(4, 25)
    B1 = frac(21, 500)

    r, s = (sqrt((15 + i * 5 * sqrt(5)) / 42) for i in [+1, -1])
    t = sqrt(frac(5, 21))

    d = {
        "zero3": [[B0]],
        "symm_rs0_roll": [[B1], [r], [s]],
        "symm_rrr": [[B1], [t]],
    }
    return S3Scheme("Ditkin 2", d, 5, _source)


def ditkin_3():
    B0 = frac(16, 175)
    B1 = frac(81, 1400)
    B2 = frac(3, 280)

    sqrt5 = sqrt(5)
    r, s = (sqrt((5 + i * sqrt5) / 18) for i in [+1, -1])
    t = sqrt(frac(1, 3))
    u, v = (sqrt((3 - i * sqrt5) / 6) for i in [+1, -1])

    d = {
        "zero3": [[B0]],
        "symm_rs0_roll": [[B1, B2], [r, u], [s, v]],
        "symm_rrr": [[B2], [t]],
    }
    return S3Scheme("Ditkin 3", d, 7, _source)


register([ditkin_1, ditkin_2, ditkin_3])
