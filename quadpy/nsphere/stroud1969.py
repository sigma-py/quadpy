# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import integrate_monomial_over_unit_nsphere
from ..helpers import untangle, pm, fsd


class Stroud1969(object):
    '''
    A.H. Stroud,
    A Fifth Degree Integration Formula for the n-Simplex,
    SIAM J. Numer. Anal., 6(1), 90–98. (9 pages),
    <https://doi.org/10.1137/0706009>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, symbolic=False):
        assert n > 2

        sqrt = numpy.vecorize(sympy.sqrt) if symbolic else numpy.sqrt
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        p_m = numpy.array([+1, -1])

        self.dim = n
        self.degree = 11

        sqrt3 = sqrt(3)

        t = sqrt(frac(1, n))
        r1, r2 = sqrt((n + 6 - p_m*4*sqrt3) / (n**2 + 12*n - 12))
        s1, s2 = sqrt((7*n - 6 + p_m*4*(n-1)*sqrt3) / (n**2 + 12*n - 12))
        u1, u2 = sqrt((n + 12 + p_m*8*sqrt3) / (n**2 + 24*n - 48))
        v1, v2 = sqrt((7*n - 12 - p_m*4*(n-2)*sqrt3) / (n**2 + 24*n - 48))

        # Solve linear equation system for x^k, k={0, 4, 6, 8, 10}, for the
        # weights (the same is done in Stroud's article).
        pts = [
            pm(n, t),
            fsd(n, (s1, 1), (r1, n-1)),
            fsd(n, (s2, 1), (r2, n-1)),
            ]
        k_range = [0, 4, 6]
        if n >= 4:
            pts.append(fsd(n, (v1, 2), (u1, n-2)))
            k_range.append(8)
        if n >= 5:
            pts.append(fsd(n, (v2, 2), (u2, n-2)))
            k_range.append(10)

        # TODO build the equation system from orthogonal polynomials
        b = [
            integrate_monomial_over_unit_nsphere([k] + (n-1)*[0])
            for k in k_range
            ]
        A = [[
            sum(p[:, 0]**k) for p in pts
            ] for k in k_range
            ]
        flt = numpy.vectorize(float)
        x = numpy.linalg.solve(flt(A), flt(b))

        data = [(x[k], pts[k]) for k in range(len(x))]

        self.points, self.weights = untangle(data)
        return
