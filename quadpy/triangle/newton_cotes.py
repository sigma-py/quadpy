# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy

from .helpers import TriangleScheme


def _newton_cotes(n, point_fun):
    """
    Construction after

    P. Silvester,
    Symmetric quadrature formulae for simplexes
    Math. Comp., 24, 95-100 (1970),
    <https://doi.org/10.1090/S0025-5718-1970-0258283-6>.
    """
    degree = n

    # points
    idx = numpy.array(
        [[i, j, n - i - j] for i in range(n + 1) for j in range(n + 1 - i)]
    )
    bary = point_fun(idx, n)

    # weights
    if n == 0:
        weights = numpy.ones(1)
        return bary, weights, degree

    def get_poly(t, m, n):
        return sympy.prod(
            [
                sympy.poly((t - point_fun(k, n)) / (point_fun(m, n) - point_fun(k, n)))
                for k in range(m)
            ]
        )

    weights = numpy.empty(len(bary))
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
    return bary, weights, degree


def NewtonCotesClosed(n):
    bary, weights, degree = _newton_cotes(n, lambda k, n: k / float(n))
    return TriangleScheme("Newton-Cotes (closed)", degree, weights, bary)


def NewtonCotesOpen(n):
    bary, weights, degree = _newton_cotes(n, lambda k, n: (k + 1) / float(n + 3))
    return TriangleScheme("Newton-Cotes (open)", degree, weights, bary)
