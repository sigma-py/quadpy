# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import orthopy
from sympy import sqrt, pi, Rational as fr

from .ditkin import Ditkin
from .hammer_stroud import HammerStroud
from .mysovskih import Mysovskih

from ..sphere import stroud as sphere_stroud
from ..helpers import untangle


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        if index == 'S3 3-1':
            self.set_data(HammerStroud('11-3'))
        elif index == 'S3 5-1':
            self.set_data(Ditkin(1))
        elif index == 'S3 5-2':
            self.set_data(Ditkin(2))
        elif index == 'S3 7-1a':
            self.set_data(HammerStroud('15-3a'))
        elif index == 'S3 7-1b':
            self.set_data(HammerStroud('15-3b'))
        elif index == 'S3 7-2':
            self.set_data(Mysovskih())
        elif index == 'S3 7-3':
            self.set_data(Ditkin(3))
        elif index == 'S3 7-4':
            # Spherical product Gauss formula.
            self.degree = 7

            # Stroud only gives decimals, sophisticated guesswork gives the
            # analytical expressions.

            # 0.9061798459, 0.5384691101
            alpha, beta = [sqrt((35 + t * 2*sqrt(70)) / 63) for t in [+1, -1]]
            rho = numpy.array([-alpha, -beta, beta, alpha])

            # 0.8611363116, 0.3399810436
            alpha, beta = [sqrt((15 + t * 2*sqrt(30)) / 35) for t in [+1, -1]]
            u = numpy.array([-alpha, -beta, beta, alpha])

            # 0.9238795325, 0.3826834324
            alpha, beta = [sqrt((2 + t * sqrt(2)) / 4) for t in [+1, -1]]
            v = numpy.array([-alpha, -beta, beta, alpha])

            # 0.1945553342, 0.1387779991
            alpha, beta = [(50 + t * sqrt(70)) / 300 for t in [+1, -1]]
            A = numpy.array([alpha, beta, beta, alpha])

            # 0.3478548451, 0.6521451549
            alpha, beta = [(18 - t * sqrt(30)) / 36 for t in [+1, -1]]
            B = numpy.array([alpha, beta, beta, alpha])

            C = numpy.full(4, pi/4)

            def outer3(a, b, c):
                '''Given 3 1-dimensional vectors a, b, c, the output is of
                shape (len(a), len(b), len(c)) and contains the values

                   out[i, j, k] = a[i] * b[j] * c[k]
                '''
                return numpy.multiply.outer(numpy.multiply.outer(a, b), c)

            vsqrt = numpy.vectorize(sqrt)
            r = outer3(rho, vsqrt(1 - u**2), vsqrt(1 - v**2))
            s = outer3(rho, vsqrt(1 - u**2), v)
            t = outer3(rho, u, 4 * [1])

            data = []
            for i in range(4):
                for j in range(4):
                    for k in range(4):
                        data.append((
                            (A[i]*B[j]*C[k]), numpy.array([[
                                r[i][j][k], s[i][j][k], t[i][j][k],
                                ]])
                            ))

            self.points, self.weights = untangle(data)
        else:
            assert index == 'S3 14-1'
            self.degree = 14

            # Get the moments corresponding to the Legendre polynomials and the
            # weight function omega(x) = x^2:
            #
            #                                    / 2/3   if k == 0,
            #    int_{-1}^{+1} |x^2| P_k(x) dx ={  8/45  if k == 2,
            #                                    \ 0     otherwise.
            #
            # In this case, the recurrence coefficients can be determined
            # analytically.
            n = 8
            alpha = numpy.full(n, fr(0))
            k = numpy.arange(n)
            beta = numpy.full(n, fr(0))
            beta[0] = fr(2, 3)
            # beta[1::2] = fr((k[1::2]+2)**2, ((2*k[1::2]+2)**2 - 1))
            for k in range(1, n, 2):
                beta[k] = fr((k+2)**2, (2*k+2)**2 - 1)
            # beta[2::2] = fr(k[2::2]**2, ((2*k[2::2]+2)**2 - 1))
            for k in range(2, n, 2):
                beta[k] = fr(k**2, (2*k+2)**2 - 1)

            # symbolic computation of the points and weights takes 4orever.
            # Keep an eye on
            # <https://math.stackexchange.com/questions/2450401/solve-small-symmetric-triadiagonal-eigenvalue-problem-symbolically>
            # for a better algorithm to be implemented in orthopy.
            flt = numpy.vectorize(float)
            alpha = flt(alpha)
            beta = flt(beta)
            points, weights = orthopy.schemes.custom(alpha, beta, mode='numpy')

            r = points[-4:]
            A = weights[-4:]

            spherical_scheme = sphere_stroud.Stroud('U3 14-1')
            v = spherical_scheme.points
            B = spherical_scheme.weights

            data = [
                (A[i]*B[j], r[i] * numpy.array([v[j]]))
                for i in range(4)
                for j in range(72)
                ]

            self.points, self.weights = untangle(data)
            self.weights *= 4.0 * numpy.pi

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
