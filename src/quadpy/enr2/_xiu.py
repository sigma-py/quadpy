import numpy as np
from sympy import Rational as frac
from sympy import cos, pi, sin, sqrt

from ..helpers import article
from ._helpers import Enr2Scheme

source = article(
    authors=["Dongbin Xiu"],
    title="Numerical integration formulas of degree two",
    journal="Applied Numerical Mathematics",
    volume="58",
    year="2008",
    pages="1515–1520",
    url="https://doi.org/10.1016/j.apnum.2007.09.004",
)


def xiu(n):
    points = []
    for k in range(n + 1):
        pt = []
        # Slight adaptation:
        # The article has points for the weight 1/sqrt(2*pi) exp(−x**2/2)
        # so divide by sqrt(2) to adapt for 1/sqrt(pi) exp(−x ** 2)
        for r in range(1, n // 2 + 1):
            alpha = (2 * r * k * pi) / (n + 1)
            pt += [cos(alpha), sin(alpha)]
        if n % 2 == 1:
            pt += [(-1) ** k / sqrt(2)]
        points.append(pt)

    points = np.array(points)
    points = np.ascontiguousarray(points.T)
    weights = np.full(n + 1, frac(1, n + 1))
    return Enr2Scheme("Xiu", n, weights, points, 2, source)
