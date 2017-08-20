# -*- coding: utf-8 -*-
#
import numpy
import orthopy

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

            # Stroud only gives decimals, sophistacated guesswork gives the
            # analytical expressions.
            plus_minus = numpy.array([+1, -1])

            rho = numpy.empty(4)
            # 0.9061798459, 0.5384691101
            rho[3], rho[2] = \
                numpy.sqrt((35.0 + plus_minus * 2*numpy.sqrt(70.0)) / 63.0)
            rho[0], rho[1] = -rho[3], -rho[2]

            u = numpy.empty(4)
            # 0.8611363116, 0.3399810436
            u[3], u[2] = \
                numpy.sqrt((15.0 + plus_minus * 2*numpy.sqrt(30.0)) / 35.0)
            u[0], u[1] = -u[3], -u[2]

            v = numpy.empty(4)
            # 0.9238795325, 0.3826834324
            v[3], v[2] = numpy.sqrt((2 + plus_minus * numpy.sqrt(2.0)) / 4.0)
            v[0], v[1] = -v[3], -v[2]

            # 0.1945553342, 0.1387779991
            A = numpy.empty(4)
            A[0], A[1] = (50.0 + plus_minus * numpy.sqrt(70.0)) / 300.0
            A[2], A[3] = A[1], A[0]

            # 0.3478548451, 0.6521451549
            B = numpy.empty(4)
            B[0], B[1] = (18.0 - plus_minus * numpy.sqrt(30.0)) / 36.0
            B[2], B[3] = B[1], B[0]

            C = numpy.full(4, numpy.pi / 4.0)

            def outer3(a, b, c):
                '''Given 3 1-dimensional vectors a, b, c, the output is of
                shape (len(a), len(b), len(c)) and contains the values

                   out[i, j, k] = a[i] * b[j] * c[k]
                '''
                return numpy.multiply.outer(numpy.multiply.outer(a, b), c)

            r = outer3(rho, numpy.sqrt(1.0 - u**2), numpy.sqrt(1.0 - v**2))
            s = outer3(rho, numpy.sqrt(1.0 - u**2), v)
            t = outer3(rho, u, numpy.ones(4))

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

            # Get the moment corresponding to the weight function omega(x) =
            # x^2. The general formula is
            #
            #                                     / 0 if k is odd,
            #    int_{-1}^{+1} |x^alpha| x^k dx ={
            #                                     \ 2/(alpha+k+1) if k is even.
            #
            n = 8
            alpha = 2.0
            k = numpy.arange(2*n+1)
            moments = (1.0 + (-1.0)**k) / (k + alpha + 1)
            gauss15 = orthopy.Gauss(n, moments)

            r = gauss15.points[-4:]
            A = gauss15.weights[-4:]

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
