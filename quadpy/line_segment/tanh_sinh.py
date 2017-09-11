# -*- coding: utf-8 -*-
#
import mpmath
from mpmath import mp


def tanh_sinh_quadrature(f, a, b, eps, f_derivatives=None):
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
    h = mp.mpf(1)
    level = 1
    while abs(error_estimate) > eps:
        # For h=1, the error estimate is too optimistic. Hence, start with
        # h=1/2 right away.
        h /= 2
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

        if f_derivatives:
            # Pretty accurate error estimation:
            #
            #   E(h) = h * (h/2/pi)**2 * sum_{-N}^{+N} F''(h*j)
            #
            # with
            #
            #   F(t) = f(g(t)) * g'(t),
            #   g(t) = tanh(pi/2 sinh(t)).
            #
            assert 1 in f_derivatives
            assert 2 in f_derivatives

            def g(t):
                return mp.tanh(mp.pi/2 * mp.sinh(t))

            def dg_dt(t):
                return mp.pi/2 * mp.cosh(t) / mp.cosh(mp.pi/2 * mp.sinh(t))**2

            def d2g_dt2(t):
                return mp.pi/2 * (
                    + mp.sinh(t)
                    - mp.pi * mp.cosh(t)**2 * mp.tanh(mp.pi/2 * mp.sinh(t))
                    ) / mp.cosh(mp.pi/2 * mp.sinh(t))**2

            def d3g_dt3(t):
                sinh_sinh = mp.sinh(mp.pi/2 * mp.sinh(t))
                cosh_sinh = mp.cosh(mp.pi/2 * mp.sinh(t))
                tanh_sinh = mp.tanh(mp.pi/2 * mp.sinh(t))
                return mp.pi/4 * mp.cosh(t) * (
                    + 2 * cosh_sinh
                    - 2 * mp.pi**2 * mp.cosh(t)**2 / cosh_sinh
                    + mp.pi**2 * mp.cosh(t)**2 * cosh_sinh
                    + mp.pi**2 * mp.cosh(t)**2 * tanh_sinh * sinh_sinh
                    - 6 * mp.pi * mp.sinh(t) * sinh_sinh
                    ) / cosh_sinh**3

            s = 0
            for jj in range(-j, +j+1):
                t = h*jj
                gt = g(t)
                g1 = dg_dt(t)
                g2 = d2g_dt2(t)
                g3 = d3g_dt3(t)
                s += (
                    + g1**3 * f_derivatives[2](gt)
                    + 3 * g1 * g2 * f_derivatives[1](gt)
                    + g3 * f(gt)
                    )
            error_estimate = s * h * (h/2/mp.pi)**2
        else:
            # more crude error estimation after Bailey
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

        print
        print('error estimate:')
        print(j, 2 - mp.fsum(weights), error_estimate)
        print

        level += 1

    exit(1)
    return value_estimates[-1], error_estimate
