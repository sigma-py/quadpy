# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._ditkin import (
    ditkin_1 as stroud_5_1,
    ditkin_2 as stroud_5_2,
    ditkin_3 as stroud_7_3,
)
from ._hammer_stroud import (
    hammer_stroud_11_3 as stroud_3_1,
    hammer_stroud_15_3a as stroud_7_1a,
    hammer_stroud_15_3b as stroud_7_1b,
)
from ._mysovskih import mysovskih as stroud_7_2

from ._helpers import BallScheme
from ..sphere import stroud as sphere_stroud
from ..helpers import untangle, book

_citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_7_4(symbolic=False):
    # spherical product gauss
    pi = sympy.pi if symbolic else numpy.pi
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    # ENH Stroud only gives decimals, sophisticated guesswork gives the analytical
    # expressions.
    pm = numpy.array([+1, -1])

    # 0.9061798459, 0.5384691101
    alpha, beta = sqrt((35 + pm * 2 * sqrt(70)) / 63)
    rho = numpy.array([-alpha, -beta, beta, alpha])

    # 0.8611363116, 0.3399810436
    alpha, beta = sqrt((15 + pm * 2 * sqrt(30)) / 35)
    u = numpy.array([-alpha, -beta, beta, alpha])

    # 0.9238795325, 0.3826834324
    alpha, beta = sqrt((2 + pm * sqrt(2)) / 4)
    v = numpy.array([-alpha, -beta, beta, alpha])

    # 0.1945553342, 0.1387779991
    alpha, beta = (50 + pm * sqrt(70)) / 300
    A = numpy.array([alpha, beta, beta, alpha])

    # 0.3478548451, 0.6521451549
    alpha, beta = (18 - pm * sqrt(30)) / 36
    B = numpy.array([alpha, beta, beta, alpha])

    C = numpy.full(4, pi / 4)

    def outer3(a, b, c):
        """Given 3 1-dimensional vectors a, b, c, the output is of
        shape (len(a), len(b), len(c)) and contains the values

           out[i, j, k] = a[i] * b[j] * c[k]
        """
        return numpy.multiply.outer(numpy.multiply.outer(a, b), c)

    r = outer3(rho, sqrt(1 - u ** 2), sqrt(1 - v ** 2))
    s = outer3(rho, sqrt(1 - u ** 2), v)
    t = outer3(rho, u, 4 * [1])

    data = [
        ((A[i] * B[j] * C[k]), numpy.array([[r[i][j][k], s[i][j][k], t[i][j][k]]]))
        for i in range(4)
        for j in range(4)
        for k in range(4)
    ]

    points, weights = untangle(data)
    return BallScheme("Stroud S3 7-4", _citation, 7, weights, points)


def stroud_14_1():
    # Get the moments corresponding to the Legendre polynomials and the weight function
    # omega(x) = x^2:
    #
    #                                    / 2/3   if k == 0,
    #    int_{-1}^{+1} |x^2| P_k(x) dx ={  8/45  if k == 2,
    #                                    \ 0     otherwise.
    #
    # In this case, the recurrence coefficients can be determined analytically.
    # ```
    # n = 8
    # alpha = numpy.full(n, fr(0))
    # k = numpy.arange(n)
    # beta = numpy.full(n, fr(0))
    # beta[0] = fr(2, 3)
    # # beta[1::2] = fr((k[1::2]+2)**2, ((2*k[1::2]+2)**2 - 1))
    # for k in range(1, n, 2):
    #     beta[k] = fr((k+2)**2, (2*k+2)**2 - 1)
    # # beta[2::2] = fr(k[2::2]**2, ((2*k[2::2]+2)**2 - 1))
    # for k in range(2, n, 2):
    #     beta[k] = fr(k**2, (2*k+2)**2 - 1)
    #
    # # symbolic computation of the points and weights takes 4orever. Keep an eye on
    # # <https://math.stackexchange.com/q/2450401/36678> for a better algorithm to be
    # # implemented in orthopy.
    # flt = numpy.vectorize(float)
    # alpha = flt(alpha)
    # beta = flt(beta)
    # points, weights = \
    #     orthopy.line.schemes.custom(alpha, beta, mode='numpy')
    #
    # r = points[-4:]
    # A = weights[-4:]
    # ```
    # TODO get symbolic expressions here
    r = numpy.array(
        [
            3.242534234038097e-01,
            6.133714327005908e-01,
            8.360311073266362e-01,
            9.681602395076261e-01,
        ]
    )
    A = numpy.array(
        [
            3.284025994586210e-02,
            9.804813271549834e-02,
            1.262636728646019e-01,
            7.618126780737085e-02,
        ]
    )

    spherical_scheme = sphere_stroud.Stroud("U3 14-1", symbolic=False)
    v = spherical_scheme.points
    B = spherical_scheme.weights

    data = [
        (A[i] * B[j], r[i] * numpy.array([v[j]])) for i in range(4) for j in range(72)
    ]

    points, weights = untangle(data)
    weights *= 4 * numpy.pi
    return BallScheme("Stroud S3 14-1", _citation, 14, weights, points)


__all__ = [
    "stroud_3_1",
    "stroud_5_1",
    "stroud_5_2",
    "stroud_7_1a",
    "stroud_7_1b",
    "stroud_7_2",
    "stroud_7_3",
    "stroud_7_4",
    "stroud_14_1",
]
