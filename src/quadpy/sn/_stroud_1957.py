import math

import numpy as np
import sympy

from ..helpers import article, untangle
from ._helpers import SnScheme

source = article(
    authors=["A.H. Stroud"],
    title="Remarks on the Disposition of Points in Numerical Integration Formulas",
    journal="Mathematical Tables and Other Aids to Computation",
    volume="11",
    number="60",
    month="oct",
    year="1957",
    pages="257-261",
    url="https://doi.org/10.2307/2001945",
)


def stroud_1957(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    sin = sympy.sin if symbolic else math.sin
    cos = sympy.cos if symbolic else math.cos
    pi = sympy.pi if symbolic else math.pi

    pts = [
        [
            [
                sqrt(frac(2, n + 2)) * cos(2 * k * i * pi / (n + 1))
                for i in range(n + 1)
            ],
            [
                sqrt(frac(2, n + 2)) * sin(2 * k * i * pi / (n + 1))
                for i in range(n + 1)
            ],
        ]
        for k in range(1, n // 2 + 1)
    ]
    if n % 2 == 1:
        sqrt3pm = np.full(n + 1, 1 / sqrt(n + 2))
        sqrt3pm[1::2] *= -1
        pts.append(sqrt3pm)
    pts = np.vstack(pts).T

    data = [(frac(1, n + 1), pts)]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Stroud 1957", n, weights, points, 2, source)
