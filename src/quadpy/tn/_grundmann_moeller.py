from math import factorial as fact

import numpy as np
from sympy import Rational as frac

from ..helpers import article, get_all_exponents, untangle
from ._helpers import TnScheme

source = article(
    authors=["A. Grundmann", "H.M. Möller"],
    title="Invariant integration formulas for the n-simplex by combinatorial methods",
    journal="SIAM J. Numer. Anal.",
    volume="15",
    year="1978",
    pages="282-290",
    url="https://doi.org/10.1137/0715019",
)


def grundmann_moeller(n: int, s: int):
    d = 2 * s + 1

    exponents = get_all_exponents(n + 1, s)

    data = [
        (
            frac(
                (-1) ** i * 2 ** (-2 * s) * (d + n - 2 * i) ** d,
                fact(i) * fact(d + n - i),
            ),
            np.array(
                [
                    [frac(2 * p + 1, d + n - 2 * i) for p in part]
                    for part in exponents[s - i]
                ]
            ),
        )
        for i in range(s + 1)
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    weights /= sum(weights)

    name = f"Grundmann-Möller(dim={n}, {s})"
    return TnScheme(name, n, weights, points, d, source)
