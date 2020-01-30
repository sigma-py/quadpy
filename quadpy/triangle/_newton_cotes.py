import math

import numpy
import sympy

from ..helpers import article
from ._helpers import TriangleScheme

citation = article(
    authors=["P. Silvester"],
    title="Symmetric quadrature formulae for simplexes",
    journal="Math. Comp.",
    volume="24",
    pages="95-100",
    year="1970",
    url="https://doi.org/10.1090/S0025-5718-1970-0258283-6",
)


def _newton_cotes(n, point_fun):
    degree = n

    # points
    idx = numpy.array(
        [[i, j, n - i - j] for i in range(n + 1) for j in range(n + 1 - i)]
    )
    points = point_fun(idx, n)

    # weights
    if n == 0:
        weights = numpy.ones(1)
        return points, weights, degree

    def get_poly(t, m, n):
        return sympy.prod(
            [
                sympy.poly((t - point_fun(k, n)) / (point_fun(m, n) - point_fun(k, n)))
                for k in range(m)
            ]
        )

    weights = numpy.empty(len(points))
    idx = 0
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            # Define the polynomial which to integrate over the
            # tetrahedron.
            t = sympy.DeferredVector("t")
            g = get_poly(t[0], i, n) * get_poly(t[1], j, n) * get_poly(t[2], k, n)
            # The integral of monomials over a tetrahedron are well-known,
            # see Silvester.
            weights[idx] = numpy.sum(
                [
                    c
                    * numpy.prod([math.factorial(l) for l in m])
                    * 2
                    / math.factorial(numpy.sum(m) + 2)
                    for m, c in zip(g.monoms(), g.coeffs())
                ]
            )
            idx += 1
    return points, weights, degree


def newton_cotes_closed(n):
    points, weights, degree = _newton_cotes(n, lambda k, n: k / float(n))
    return TriangleScheme(
        f"Newton-Cotes (closed, {n})", weights, points, degree, citation
    )


def newton_cotes_open(n):
    points, weights, degree = _newton_cotes(n, lambda k, n: (k + 1) / float(n + 3))
    if n == 0:
        degree = 1
    return TriangleScheme(
        f"Newton-Cotes (open, {n})", weights, points, degree, citation
    )
