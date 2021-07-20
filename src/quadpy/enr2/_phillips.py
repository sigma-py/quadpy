import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, untangle
from ._helpers import Enr2Scheme

source = article(
    authors=["G.M. Phillips"],
    title="A survey of one-dimensional and multidimensional numerical integration",
    journal="Computer Physics Communications",
    volume="20",
    number="1",
    year="1980",
    pages="17-27",
    url="https://doi.org/10.1016/0010-4655%2880%2990102-2",
)


def phillips(n):
    lmbda = sqrt(frac(3, 2))
    w1 = frac(n ** 2 - 7 * n + 18, 18)
    w2 = frac(4 - n, 18)
    w3 = frac(1, 36)

    data = [
        (w1, np.full((1, n), 0)),
        (w2, fsd(n, (lmbda, 1))),
        (w3, fsd(n, (lmbda, 2))),
    ]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Phillips", n, weights, points, 5, source)
