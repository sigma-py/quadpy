import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, untangle, z
from ._helpers import SnScheme

source = article(
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
    A = frac(1, 8)
    B = frac(5 - n, 4)
    C = frac((6 - n) * (1 - n ** 2) + 36, 4 * (n + 3))
    D = frac(81, (n + 3) * (n + 6) ** 2)
    E = frac(45 * n ** 2 + 324 * n + 216, n ** 2 + 12 * n + 36) - frac(
        n * (n ** 2 - 12 * n + 65), 6
    )

    r = sqrt(frac(3, n + 6))
    data = [
        (A, fsd(n, (r, 3))),
        (B, fsd(n, (r, 2))),
        (C, fsd(n, (r, 1))),
        (D, fsd(n, (1, 1))),
        (E, z(n)),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)

    weights /= frac(27 * (n + 2) * (n + 4), (n + 6) ** 2)
    return SnScheme("Dobrodeev 1970", n, weights, points, 7, source)
