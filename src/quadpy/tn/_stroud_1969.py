import ndim
import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, rd, untangle
from ._helpers import TnScheme

source = article(
    authors=["A.H. Stroud"],
    title="A Fifth Degree Integration Formula for the n-Simplex",
    journal="SIAM J. Numer. Anal.",
    volume="6",
    number="1",
    pages="90–98",
    url="https://doi.org/10.1137/0706009",
)


def stroud_1969(n):
    assert n >= 3

    degree = 5

    sqrt15 = sqrt(15)

    t = frac(1, n + 1)
    r1, r2 = ((n + 4 - pm * sqrt15) / (n ** 2 + 8 * n + 1) for pm in [+1, -1])
    s1, s2 = ((4 * n + 1 + pm * n * sqrt15) / (n ** 2 + 8 * n + 1) for pm in [+1, -1])
    u1, u2 = ((n + 7 + pm * 2 * sqrt15) / (n ** 2 + 14 * n - 11) for pm in [+1, -1])
    v1, v2 = (
        (4 * n - 2 - pm * (n - 1) * sqrt15) / (n ** 2 + 14 * n - 11) for pm in [+1, -1]
    )

    # Solve linear equation system for x^k, k={0, 2, 3, 4, 5}, for the
    # weights (the same is done in Stroud's article).
    pts = [
        np.full((1, n + 1), t),
        rd(n + 1, [(r1, n), (s1, 1)]),
        rd(n + 1, [(r2, n), (s2, 1)]),
        rd(n + 1, [(u1, n - 1), (v1, 2)]),
    ]
    k_range = [0, 2, 3, 4]

    if n > 3:
        pts.append(rd(n + 1, [(u2, n - 1), (v2, 2)]))
        k_range.append(5)

    b0 = ndim.nsimplex.integrate_monomial(n * [0], symbolic=True)
    b = [
        ndim.nsimplex.integrate_monomial(np.array([k] + (n - 1) * [0]), symbolic=True)
        / b0
        for k in k_range
    ]

    A = [[sum(p[:, 0] ** k) for p in pts] for k in k_range]

    flt = np.vectorize(float)
    w = np.linalg.solve(flt(A), flt(b))

    data = [(w[i], pts[i]) for i in range(len(w))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return TnScheme("Stroud 1969", n, weights, points, degree, source)
