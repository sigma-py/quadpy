# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import integrate_monomial_over_unit_simplex
from ..helpers import untangle, rd


class Stroud1969(object):
    '''
    A.H. Stroud,
    A Fifth Degree Integration Formula for the n-Simplex,
    SIAM J. Numer. Anal., 6(1), 90–98. (9 pages),
    <https://doi.org/10.1137/0706009>.
    '''

    def __init__(self, n, symbolic=False):
        assert n >= 3

        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        self.dim = n
        self.degree = 5

        sqrt15 = sqrt(15)

        t = frac(1, n+1)
        r1, r2 = [(n + 4 - pm * sqrt15) / (n**2 + 8*n + 1) for pm in [+1, -1]]
        s1, s2 = [
            (4*n + 1 + pm * n * sqrt15) / (n**2 + 8*n + 1) for pm in [+1, -1]
            ]
        u1, u2 = [
            (n + 7 + pm * 2 * sqrt15) / (n**2 + 14*n - 11) for pm in [+1, -1]
            ]
        v1, v2 = [
            (4*n - 2 - pm * (n-1) * sqrt15) / (n**2 + 14*n - 11)
            for pm in [+1, -1]
            ]

        # Solve linear equation system for x^k, k={0, 2, 3, 4, 5}, for the
        # weights (the same is done in Stroud's article).
        pts = [
            numpy.full((1, n+1), t),
            rd(n+1, [(r1, n), (s1, 1)]),
            rd(n+1, [(r2, n), (s2, 1)]),
            rd(n+1, [(u1, n-1), (v1, 2)]),
            ]
        k_range = [0, 2, 3, 4]

        if n > 3:
            pts.append(
                rd(n+1, [(u2, n-1), (v2, 2)])
                )
            k_range.append(5)

        b0 = integrate_monomial_over_unit_simplex(n*[0], symbolic=symbolic)
        b = [
            integrate_monomial_over_unit_simplex(
                numpy.array([k] + (n-1)*[0]), symbolic=symbolic
                ) / b0
            for k in k_range
            ]

        A = [[
            sum(p[:, 0]**k) for p in pts
            ] for k in k_range
            ]

        flt = numpy.vectorize(float)
        x = numpy.linalg.solve(flt(A), flt(b))

        data = [(x[i], pts[i]) for i in range(len(x))]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
