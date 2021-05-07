import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, comb, fsd, untangle, z
from ._helpers import CnScheme

_source = article(
    authors=["G.M. Phillips"],
    title="Numerical integration over an N-dimensional rectangular region",
    journal="Comput J",
    year="1967",
    volume="10",
    number="3",
    pages="297-299",
    url="https://doi.org/10.1093/comjnl/10.3.297",
)


def phillips(n):
    if n == 2:
        p1 = 1
        p2 = frac(14, 3)
        q = frac(5, 3)
    elif n == 3:
        p1 = 1
        p2 = frac(14, 5)
        q = frac(5, 2)
        r = 1
    elif n == 4:
        p1 = 1
        p2 = frac(112, 11)
        q = 5
        r = 2
    else:
        assert n >= 5
        p1 = 1
        En = frac(25 * n ** 2 - 165 * n + 302, 972)
        p2 = 1 / (frac(3, 5) - frac(1, 35 * En))
        q = frac(5, 3)
        r = frac(5, 3)

    gamma = frac((n - 1) * (19 - 5 * n), 270)
    delta = frac((n - 1) * (n - 2), 108)

    a1 = frac(23 - 5 * n, 180) - gamma * q / 2
    a2 = frac(35 * n ** 2 - 231 * n + 466, 3780)
    beta1 = (a1 - a2 * p2) / (p1 - p2)
    beta2 = (a1 - a2 * p1) / (p2 - p1)

    lambda1 = 1 / sqrt(p1)
    lambda2 = 1 / sqrt(p2)
    mu = 1 / sqrt(q)

    b1 = beta1 / lambda1 ** 6
    b2 = beta2 / lambda2 ** 6

    c = gamma / (2 * (n - 1) * mu ** 6)

    a = 1 - 2 * n * (b1 + b2) - 4 * comb(n, 2) * c

    if n > 2:
        nu = 1 / sqrt(r)
        d = delta / (4 * comb(n - 1, 2) * nu ** 6)
        a -= 8 * comb(n, 3) * d

    data = [
        (a, z(n)),
        (b1, fsd(n, (lambda1, 1))),
        (b2, fsd(n, (lambda2, 1))),
        (c, fsd(n, (mu, 2))),
    ]

    if n > 2:
        data.append((d, fsd(n, (nu, 3))))

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return CnScheme("Phillips", n, weights, points, 7, _source, 1.521e-13)
