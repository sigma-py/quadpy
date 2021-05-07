import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, get_nsimplex_points
from ._helpers import UnScheme

source = article(
    authors=["I.P. Mysovskikh"],
    title="The approximation of multiple integrals by using interpolatory cubature formulae",
    journal="Proceedings of a Symposium on Quantitative Approximation Held in Bonn, West Germany, August 20–24, 1979",
    year="1980",
    pages="217-243",
    url="https://doi.org/10.1016/B978-0-12-213650-4.50025-8",
)


def mysovskikh_1(n):
    points = get_nsimplex_points(n, sqrt, frac)
    weights = np.full(n + 1, frac(1, n + 1))
    return UnScheme("Mysovskikh 1", n, weights, points, 2, source)


def mysovskikh_2(n):
    A = frac((7 - n) * n, 2 * (n + 1) ** 2 * (n + 2))
    B = frac(2 * (n - 1) ** 2, n * (n + 1) ** 2 * (n + 2))

    a = get_nsimplex_points(n, sqrt, frac)
    b = np.array(
        [
            sqrt(frac(n, 2 * (n - 1))) * (a[k] + a[l])
            for k in range(n + 1)
            for l in range(k)
        ]
    )

    points = np.concatenate([a, -a, b, -b])
    weights = np.concatenate(
        [
            np.full(len(a), A),
            np.full(len(a), A),
            np.full(len(b), B),
            np.full(len(b), B),
        ]
    )
    return UnScheme("Mysovskikh 2", n, weights, points, 5, source)
