import math

import sympy

from ..helpers import article
from ._helpers import C3Scheme, register

_source = article(
    authors=["Preston C. Hammer", "A. Wayne Wymore"],
    title="Numerical evaluation of multiple integrals. I",
    journal="Math. Comp.",
    volume="11",
    year="1957",
    pages="59-67",
    url="https://doi.org/10.1090/S0025-5718-1957-0087220-6",
)


# hw(27 / 20) has all the desirable properties (positive weights, points inside the
# domain) while minimizing the ratio max(weights) / min(weights)
def hammer_wymore(lmbda=sympy.Rational(27, 20)):
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

    b = sqrt(c1 ** 2 - 4 * c0)
    u3 = (-c1 - b) / 2
    u4 = (-c1 + b) / 2

    a3 = (p0 * u4 - p1) / 8 / (u4 - u3)
    a4 = (p0 * u3 - p1) / 8 / (u3 - u4)

    x1 = sqrt(u1)
    x2 = sqrt(u2)
    x3 = sqrt(u3)
    x4 = sqrt(u4)

    d = {
        "symm_r00": [[a1 / 8], [x1]],
        "symm_rr0": [[a2 / 8], [x2]],
        "symm_rrr": [[a3 / 8, a4 / 8], [x3, x4]],
    }
    return C3Scheme(f"Hammer-Wymore (lambda = {lmbda})", d, 7, _source)


register([hammer_wymore])
