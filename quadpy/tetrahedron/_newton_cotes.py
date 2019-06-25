# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy

from ._helpers import TetrahedronScheme
from ..helpers import article

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
        [
            [i, j, k, n - i - j - k]
            for i in range(n + 1)
            for j in range(n + 1 - i)
            for k in range(n + 1 - i - j)
        ]
    )
    points = point_fun(idx, n)

    # weights
    if n == 0:
        weights = numpy.ones(1)
        return points, points, weights, degree

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
            for k in range(n + 1 - i - j):
                L = n - i - j - k
                # Compute weight.
                # Define the polynomial which to integrate over the
                # tetrahedron.
                t = sympy.DeferredVector("t")
                g = (
                    get_poly(t[0], i, n)
                    * get_poly(t[1], j, n)
                    * get_poly(t[2], k, n)
                    * get_poly(t[3], L, n)
                )
                # The integral of monomials over a tetrahedron are well-known,
                # see Silvester.
                weights[idx] = numpy.sum(
                    [
                        c
                        * numpy.prod([math.factorial(k) for k in m])
                        * 6.0
                        / math.factorial(numpy.sum(m) + 3)
                        for m, c in zip(g.monoms(), g.coeffs())
                    ]
                )
                idx += 1
    return weights, points, degree, citation


def newton_cotes_closed(n):
    return TetrahedronScheme(
        "Newton-Cotes (closed, {})".format(n),
        *_newton_cotes(n, lambda k, n: k / float(n))
    )


def newton_cotes_open(n):
    scheme = TetrahedronScheme(
        "Newton-Cotes (open, {})".format(n),
        *_newton_cotes(n, lambda k, n: (k + 1) / float(n + 4))
    )
    if n == 0:
        scheme.degree = 1
    return scheme
