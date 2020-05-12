import math

import sympy

from ..helpers import article, untangle
from ._helpers import HexahedronScheme, fs_r00, fs_rr0, pm_rrr

_citation = article(
    authors=["Preston C. Hammer", "A. Wayne Wymore"],
    title="Numerical evaluation of multiple integrals. I",
    journal="Math. Comp.",
    volume="11",
    year="1957",
    pages="59-67",
    url="https://doi.org/10.1090/S0025-5718-1957-0087220-6",
)


def hammer_wymore(lmbda=1):
    symbolic = not isinstance(lmbda, float)
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else math.sqrt

    I0 = 8
    I2 = frac(8, 3)
    I4 = frac(8, 5)
    I22 = frac(8, 9)
    I6 = frac(8, 7)
    I42 = frac(8, 15)
    I222 = frac(8, 27)

    u1 = (I6 - 2 * I42 + I222 + lmbda * (I42 - I222)) / (I4 - I22)
    u2 = u1 / lmbda
    a1 = (I6 - 2 * I42 + I222) / 2 / u1 ** 3
    a2 = (I42 - I222) / 4 / u2 ** 3

    p0 = I0 - 6 * a1 - 12 * a2
    p1 = I2 - 2 * a1 * u1 - 8 * a2 * u2
    p2 = I22 - 4 * a2 * u2 ** 2
    p3 = I222

    c0 = (p3 * p1 - p2 ** 2) / (p0 * p2 - p1 ** 2)
    c1 = (p3 * p0 - p1 * p2) / (p1 ** 2 - p2 * p0)

    u3 = (-c1 - sqrt(c1 ** 2 - 4 * c0)) / 2
    u4 = (-c1 + sqrt(c1 ** 2 - 4 * c0)) / 2

    a3 = (p0 * u4 - p1) / 8 / (u4 - u3)
    a4 = (p0 * u3 - p1) / 8 / (u3 - u4)

    x1 = sqrt(u1)
    x2 = sqrt(u2)
    x3 = sqrt(u3)
    x4 = sqrt(u4)

    data = [(a1, fs_r00(x1)), (a2, fs_rr0(x2)), (a3, pm_rrr(x3)), (a4, pm_rrr(x4))]

    points, weights = untangle(data)
    return HexahedronScheme(
        f"Hammer-Wymore (lambda = {lmbda})", weights, points, 7, _citation
    )
