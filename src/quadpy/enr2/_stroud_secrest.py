import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, pm, untangle
from ._helpers import Enr2Scheme

source = article(
    authors=["A.H. Stroud", "D. Secrest"],
    title="Approximate integration formulas for certain spherically symmetric regions",
    journal="Math. Comp.",
    volume="17",
    year="1963",
    pages="105-135",
    url="https://doi.org/10.1090/S0025-5718-1963-0161473-0",
)


def stroud_secrest_1(n):
    # TODO check which is more appropriate
    # print(_nsimplex(n))
    # print()
    # print(get_nsimplex_points(n))
    data = [(frac(1, n + 1), sqrt(frac(1, 2)) * _nsimplex(n))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud-Secrest I", n, weights, points, 2, source)


def stroud_secrest_2(n):
    nu = sqrt(frac(n, 2))
    data = [(frac(1, 2 * n), fsd(n, (nu, 1)))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud-Secrest II", n, weights, points, 3, source)


def stroud_secrest_3(n):
    nu = sqrt(frac(1, 2))
    data = [(frac(1, 2 ** n), pm(n * [nu]))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud-Secrest III", n, weights, points, 3, source)


def stroud_secrest_4(n):
    nu = sqrt(frac(n + 2, 2))
    xi = sqrt(frac(n + 2, 4))
    A = frac(2, n + 2)
    B = frac(4 - n, 2 * (n + 2) ** 2)
    C = frac(1, (n + 2) ** 2)

    data = [(A, np.full((1, n), 0)), (B, fsd(n, (nu, 1))), (C, fsd(n, (xi, 2)))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud-Secrest IV", n, weights, points, 5, source)


def _nsimplex(n):
    # construct the regular n-simplex points with 0 center
    return np.array(
        [
            [-sqrt(frac(n + 1, (n + 1 - k) * (n - k))) for k in range(i)]
            + [sqrt(frac((n + 1) * (n - i), n + 1 - i))]
            + (n - i - 1) * [0]
            for i in range(n)
        ]
        + [[-sqrt(frac(n + 1, (n + 1 - i) * (n - i))) for i in range(n)]]
    )
