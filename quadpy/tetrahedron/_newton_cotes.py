import math

import numpy
import sympy

from ..helpers import article, prod, get_all_exponents
from ._helpers import TetrahedronScheme

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
    idxs = numpy.array(get_all_exponents(3 + 1, n)[-1])
    points = point_fun(idxs, n)

    # weights
    if n == 0:
        weights = numpy.ones(1)
        return weights, points, weights, degree

    def get_poly(t, m, n):
        return sympy.prod(
            [
                sympy.poly((t - point_fun(k, n)) / (point_fun(m, n) - point_fun(k, n)))
                for k in range(m)
            ]
        )

    weights = numpy.empty(len(points))
    kk = 0
    d = 3
    for idx in idxs:
        # Compute weight.
        # Define the polynomial which to integrate over the tetrahedron.
        t = sympy.DeferredVector("t")
        g = prod(get_poly(t[k], i, n) for k, i in enumerate(idx))
        # The integral of monomials over a tetrahedron are well-known, see Silvester.
        weights[kk] = numpy.sum(
            [
                c
                * prod([math.factorial(k) for k in m])
                * math.factorial(d)
                / math.factorial(numpy.sum(m) + d)
                for m, c in zip(g.monoms(), g.coeffs())
            ]
        )
        kk += 1
    return weights, points, degree, citation


def newton_cotes_closed(n):
    return TetrahedronScheme(
        f"Newton-Cotes (closed, {n})", *_newton_cotes(n, lambda k, n: k / float(n))
    )


def newton_cotes_open(n):
    n = 3
    scheme = TetrahedronScheme(
        f"Newton-Cotes (open, {n})",
        *_newton_cotes(n, lambda k, n: (k + 1) / float(n + 4)),
    )
    print(scheme.points)
    print(scheme.weights)
    exit(1)
    if n == 0:
        scheme.degree = 1
    return scheme
