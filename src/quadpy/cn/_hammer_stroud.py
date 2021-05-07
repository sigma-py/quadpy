import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, expand_symmetries, fsd, untangle, z
from ._helpers import CnScheme

_source = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    volume="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)


def hammer_stroud_1n(n):
    d = {"a0": [[frac(1, 2 * n)], [sqrt(frac(n, 3))]]}
    points, weights = expand_symmetries(d, n)
    return CnScheme("Hammer-Stroud 1n", n, weights, points, 3, _source, 5.863e-14)


def hammer_stroud_2n(n):
    r = sqrt(frac(3, 5))
    data = [
        (frac(25 * n ** 2 - 115 * n + 162, 162), z(n)),
        (frac(70 - 25 * n, 162), fsd(n, (r, 1))),
        (frac(25, 324), fsd(n, (r, 2))),
    ]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return CnScheme("Hammer-Stroud 2n", n, weights, points, 5, _source, 3.820e-14)
