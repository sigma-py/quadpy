import math

import numpy as np
import sympy

from ..helpers import article, combine, fsd, pm, untangle, z
from ._helpers import SnScheme

source = article(
    authors=["A.H. Stroud"],
    title="Some Fifth Degree Integration Formulas for Symmetric Regions",
    journal="Mathematics of Computation",
    volume="20",
    number="93",
    month="jan",
    year="1966",
    pages="90-97",
    url="https://doi.org/10.1090/S0025-5718-1966-0191094-8",
)


def stroud_1966_a(n, symbolic=False):
    sqrt = sympy.sqrt if symbolic else math.sqrt

    a = sqrt(2 * (n + 4))

    r2 = (n + 4 - a) / (n + 4)
    s2 = (n * (n + 4) + 2 * a) / (n ** 2 + 2 * n - 4) / (n + 4)

    B1 = 1 / r2 ** 2 / (n + 2) / (n + 4)
    B2 = 1 / s2 ** 2 / 2 ** n / (n + 2) / (n + 4)

    r = sqrt(r2)
    s = sqrt(s2)

    data = [(B1, fsd(n, (r, 1))), (B2, pm(n * [s]))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Stroud 1966a", n, weights, points, 5, source)


def stroud_1966_b(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt

    alpha = 0
    s = sqrt(frac(n + alpha + 2, (n + 2) * (n + alpha + 4)))
    data = []
    B0 = 1
    for k in range(1, n + 1):
        B = frac(
            2 ** (k - n) * (n + 2) * (n + alpha) * (n + alpha + 4),
            n * (k + 1) * (k + 2) * (n + alpha + 2) ** 2,
        )
        B0 -= 2 ** (n - k + 1) * B
        r = sqrt(frac((k + 2) * (n + alpha + 2), (n + 2) * (n + alpha + 4)))
        v = np.concatenate(
            [
                np.zeros((2 ** (n - k + 1), k - 1), dtype=int),
                pm(np.array([r] + (n - k) * [s])),
            ],
            axis=-1,
        )
        data.append((B, v))
    data.append((B0, z(n)))

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Stroud 1966b", n, weights, points, 5, source)


def stroud_1966_c(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt

    a = sqrt(2 * (n + 2))
    r = sqrt((n + 2 + (n - 1) * a) / (n * (n + 4)))
    s = sqrt((n + 2 - a) / (n * (n + 4)))

    B0 = frac(4, (n + 2) ** 2)
    B1 = frac(n + 4, 2 ** n * (n + 2) ** 2)

    data = [(B0, z(n)), (B1, combine(((+r, -r), 1), ((+s, -s), (n - 1))))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Stroud 1966c", n, weights, points, 5, source)


def stroud_1966_d(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt

    a = 2 * sqrt(n + 4)
    b = sqrt(2 * (n + 1) * (n + 2) * (n + 4))

    r = sqrt((n * (n + 4) + a + (n - 1) * b) / (n * (n + 2) * (n + 4)))
    s = sqrt((n * (n + 4) + a - b) / (n * (n + 2) * (n + 4)))
    t = sqrt((n + 4 - a) / ((n + 2) * (n + 4)))

    B = frac(1, 2 ** n * (n + 1))

    # The data is given symbolically, and for large n, those are thousands of points and
    # weights. Converting them to float takes a long time. A better approach would be to
    # be convert r, s, t first, and assemble the data afterwards.
    data = [(B, combine(((+r, -r), 1), ((+s, -s), (n - 1)))), (B, pm(n * [t]))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Stroud 1966d", n, weights, points, 5, source, 1.019e-14)
