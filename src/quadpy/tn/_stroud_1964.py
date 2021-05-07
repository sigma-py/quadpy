import numpy as np
from mpmath import mp

from ..helpers import article, rd, untangle
from ._helpers import TnScheme

source = article(
    authors=["A.H. Stroud"],
    title="Approximate Integration Formulas of Degree 3 for Simplexes",
    journal="Mathematics of Computation",
    volume="18",
    number="88",
    month="oct",
    year="1964",
    pages="590-597",
    url="https://doi.org/10.2307/2002945",
)


def _stroud_1964(variant_a, n):
    degree = 3

    if n == 2:
        # The roots sum up to 1
        r, s, t = mp.polyroots([1, -1, 0.25, -1.0 / 60.0])
        data = [(1.0 / (n * (n + 1)), rd(n + 1, [(r, 1), (s, 1), (t, 1)]))]
    else:
        assert n > 2

        # Stroud's book only gives numerical values for certain n; the article explains
        # it in more detail, namely: r is a root of a polynomial of degree 3.
        rts = np.sort(np.roots([n + 1, -3, 3 / (n + 2), -1 / ((n + 2) * (n + 3))]))

        # all roots are real-valued
        if n > 8:
            assert not variant_a, "Choose variant b for n >= 9."

        r = rts[0] if variant_a else rts[1]

        # s and t are zeros of a polynomial of degree 2
        s, t = np.sort(
            np.roots(
                [
                    1,
                    -(1 - (n - 1) * r),
                    n / (2 * (n + 2)) - (n - 1) * r + n * (n - 1) / 2 * r ** 2,
                ]
            )
        )

        data = [(1 / (n * (n + 1)), rd(n + 1, [(r, n - 1), (s, 1), (t, 1)]))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)

    variant = "a" if variant_a else "b"
    name = f"Stroud 1964{variant}"
    return TnScheme(name, n, weights, points, degree, source)


def stroud_1964a(n):
    return _stroud_1964(True, n)


def stroud_1964b(n):
    return _stroud_1964(False, n)
