# -*- coding: utf-8 -*-
#
import mpmath
from mpmath import mp


def tanh_sinh_quadrature(f, a, b, eps):
    '''Integrate a function `f` between `a` and `b` with accuracy `eps`.


    Mori, Masatake
    Discovery of the double exponential transformation and its developments,
    Publications of the Research Institute for Mathematical Sciences,
    41 (4): 897–935, ISSN 0034-5318,
    doi:10.2977/prims/1145474600,
    <http://www.kurims.kyoto-u.ac.jp/~okamoto/paper/Publ_RIMS_DE/41-4-38.pdf>.
    '''
    # David H. Bailey, Karthik Jeyabalan, and Xiaoye S. Li,
    # Error function quadrature,
    # Experiment. Math., Volume 14, Issue 3 (2005), 317-329,
    # <https://projecteuclid.org/euclid.em/1128371757>.
    #
    # David H. Bailey,
    # Tanh-Sinh High-Precision Quadrature,
    # 2006,
    # <http://www.davidhbailey.com/dhbpapers/dhb-tanh-sinh.pdf>.
    assert a < b

    # ts = mpmath.calculus.quadrature.TanhSinh(mpmath.mp)
    # print(ts)
    # degree = 3
    # points, weights = zip(*ts.calc_nodes(degree, 10))
    # import numpy
    # # import matplotlib.pyplot as plt
    # # x = numpy.linspace(-10, +10, 10000)
    # # # plt.plot(x, numpy.tanh(numpy.pi/2 * numpy.sinh(x)))
    # # plt.semilogy(
    # #     x,
    # #     numpy.pi/2 * numpy.cosh(x) / numpy.cosh(numpy.pi/2 * numpy.sinh(x))**2
    # #     )
    # # plt.show()
    # # exit(1)
    # n = 4
    # h = 1.0 / 2**3
    # j = 0
    # weights = [h * mp.pi/2]
    # u2 = []
    # cosh_u2 = []
    # while weights[-1] > 1.0e-40:
    #     j += 1
    #     u1 = h * mp.pi/2 * mp.cosh(h*j)
    #     u2.append(mp.pi/2 * mp.sinh(h*j))
    #     cosh_u2.append(mp.cosh(u2[-1]))
    #     weights.append(u1 / cosh_u2[-1]**2)
    # # mirror
    # weights += weights[1:]
    # print(len(weights))
    # print(2 - mp.fsum(weights))
    # exit(1)

    # n = 3
    # h = 1.0 / 2**2
    # hj = numpy.linspace(-n, n, int(2*n/h)+1)
    # print(hj)
    # weights = (
    #     h * numpy.pi/2 * numpy.cosh(hj)
    #     / numpy.cosh(numpy.pi/2 * numpy.sinh(hj))**2
    #     )
    # print(weights)
    # print(weights[0])
    # import math
    # print(2 - math.fsum(weights))

    # # exit(1)

    num_digits = int(-mp.log10(eps) + 1)
    mpmath.mp.dps = num_digits

    value_estimates = []
    error_estimate = mp.mpf(1)
    h = 1
    level = 1
    while error_estimate > eps:
        # Compute all weights until w < eps**2
        weights = [h * mp.pi/2]
        j = 0
        u2 = []
        cosh_u2 = []
        while weights[-1] > eps**2:
            j += 1
            u1 = h * mp.pi/2 * mp.cosh(h*j)
            u2.append(mp.pi/2 * mp.sinh(h*j))
            cosh_u2.append(mp.cosh(u2[-1]))
            weights.append(u1 / cosh_u2[-1]**2)
        # mirror
        weights += weights[1:]

        print(h, h*j, len(weights), 2 - mp.fsum(weights))

        # x = [mp.tanh(v) for v in pi2_sinh_hj]
        # y = 1 - x
        y = [1 / (mp.exp(v) * c) for v, c in zip(u2, cosh_u2)]
        # mirror at the origin
        y += [-yy for yy in y[1:]]

        # perform the integration
        summands = [f(1 - yy) * weight for yy, weight in zip(y, weights)]
        value_estimates.append(mp.fsum(summands))

        # print(value_estimates[-1])
        # print

        # error estimation after Bailey
        if level <= 2:
            error_estimate = 1.0
        elif value_estimates[0] == value_estimates[-1]:
            error_estimate = 0.0
        else:
            d1 = mp.log10(abs(value_estimates[0] - value_estimates[-1]))
            d2 = mp.log10(abs(value_estimates[0] - value_estimates[-2]))
            d3 = mp.log10(eps * max([abs(x) for x in summands]))
            d4 = mp.log10(max(summands[0], summands[-1]))
            d = int(max(d1**2 / d2, 2*d1, d3, d4))
            # print(d1)
            # print(d2)
            # print(d3)
            # print(d4)
            # print(d)
            # print
            error_estimate = 10**d

        h /= 2.0
        level += 1

    exit(1)
    return value_estimates[-1], error_estimate
