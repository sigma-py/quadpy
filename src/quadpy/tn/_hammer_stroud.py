import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, rd, untangle
from ._helpers import TnScheme

source = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Integration Over Simplexes",
    journal="Mathematical Tables and Other Aids to Computation",
    volume="10",
    number="55",
    month="jul",
    year="1956",
    pages="137-139",
    url="https://doi.org/10.2307/2002484",
)


def hammer_stroud_1a(n: int):
    degree = 2
    r = (n + 2 - sqrt(n + 2)) / (n + 1) / (n + 2)
    s = (n + 2 + n * sqrt(n + 2)) / (n + 1) / (n + 2)
    data = [(frac(1, n + 1), rd(n + 1, [(r, n), (s, 1)]))]

    points, weights = untangle(data)
    return TnScheme("Hammer-Stround 1a", n, weights, points, degree, source)


def hammer_stroud_1b(n: int):
    degree = 2
    r = (n + 2 + sqrt(n + 2)) / (n + 1) / (n + 2)
    s = (n + 2 - n * sqrt(n + 2)) / (n + 1) / (n + 2)
    data = [(frac(1, n + 1), rd(n + 1, [(r, n), (s, 1)]))]

    points, weights = untangle(data)
    return TnScheme("Hammer-Stround 1b", n, weights, points, degree, source)


def hammer_stroud_2(n: int):
    degree = 3

    B = -frac((n + 1) ** 2, 4 * (n + 2))
    C = frac((n + 3) ** 2, 4 * (n + 1) * (n + 2))

    r = frac(1, n + 1)
    s = frac(1, n + 3)
    t = frac(3, n + 3)

    data = [(B, [(n + 1) * [r]]), (C, rd(n + 1, [(t, 1), (s, n)]))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return TnScheme("Hammer-Stround 2", n, weights, points, degree, source)
