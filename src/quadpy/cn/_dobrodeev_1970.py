import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, untangle, z
from ._helpers import CnScheme

_source = article(
    authors=["L.N. Dobrodeev"],
    title="Cubature formulas of the seventh order of accuracy for a hypersphere and a hypercube",
    journal="USSR Computational Mathematics and Mathematical Physics",
    volume="10",
    number="1",
    year="1970",
    pages="252–253",
    url="https://doi.org/10.1016/0041-5553%2870%2990084-4",
)


def dobrodeev_1970(n):
    assert n >= 3, f"Only works for n >= 3, not n = {n}"

    A = frac(1, 8)
    B = frac(19 - 5 * n, 20)
    alpha = 35 * n * (5 * n - 33)
    C = frac((alpha + 2114) ** 3, 700 * (alpha + 1790) * (alpha + 2600))
    D = frac(729, 1750) * frac(alpha + 2114, alpha + 2600)
    E = frac(n * (n - 1) * (n - frac(47, 10)), 3) - 2 * n * (C + D) + frac(729, 125)

    a = sqrt(frac(3, 5))
    b = a
    c = sqrt(frac(3, 5) * frac(alpha + 1790, alpha + 2114))
    d = 1
    data = [
        (A, fsd(n, (a, 3))),
        (B, fsd(n, (b, 2))),
        (C, fsd(n, (c, 1))),
        (D, fsd(n, (d, 1))),
        (E, z(n)),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    weights *= frac(125, 729)

    # d = {
    #     "0": [[w0]],
    #     "a0": [[w1], [lmbdas[0]]],
    #     "aa0": [[w1], [lmbdas[0]]],
    #     "aaa0": [[w1], [lmbdas[0]]],
    #     "a": [[w], [sqrt(delta2)]]
    # }
    return CnScheme("Dobrodeev 1970", n, weights, points, 7, _source, 8.100e-13)
