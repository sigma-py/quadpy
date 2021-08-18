import ndim
import numpy as np
from numpy import sqrt

from ..helpers import article, fsd, pm, untangle
from ._helpers import UnScheme

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
    assert n > 2
    degree = 11

    sqrt3 = sqrt(3)

    t = sqrt(1 / n)
    r1, r2 = (sqrt((n + 6 - i * 4 * sqrt3) / (n ** 2 + 12 * n - 12)) for i in [+1, -1])
    s1, s2 = (
        sqrt((7 * n - 6 + i * 4 * (n - 1) * sqrt3) / (n ** 2 + 12 * n - 12))
        for i in [+1, -1]
    )
    u1, u2 = (sqrt((n + 12 + i * 8 * sqrt3) / (n ** 2 + 24 * n - 48)) for i in [+1, -1])
    v1, v2 = (
        sqrt((7 * n - 12 - i * 4 * (n - 2) * sqrt3) / (n ** 2 + 24 * n - 48))
        for i in [+1, -1]
    )

    # Solve linear equation system for x^k, k={0, 4, 6, 8, 10}, for the weights (the
    # same is done in Stroud's article).
    pts = [pm(n * [t]), fsd(n, (s1, 1), (r1, n - 1)), fsd(n, (s2, 1), (r2, n - 1))]
    k_range = [0, 4, 6]
    if n >= 4:
        pts.append(fsd(n, (v1, 2), (u1, n - 2)))
        k_range.append(8)
    if n >= 5:
        pts.append(fsd(n, (v2, 2), (u2, n - 2)))
        k_range.append(10)
    # TODO build the equation system from orthogonal polynomials
    b = [
        ndim.nsphere.integrate_monomial([k] + (n - 1) * [0], symbolic=False)
        for k in k_range
    ]
    A = [[sum(p[:, 0] ** k) for p in pts] for k in k_range]
    w = np.linalg.solve(A, b)

    data = [(w[k], pts[k]) for k in range(len(w))]

    points, weights = untangle(data)
    weights /= ndim.nsphere.volume(n)
    return UnScheme("Stroud 1969", n, weights, points, degree, source, 4.270e-14)
