# -*- coding: utf-8 -*-
#
import mpmath
from mpmath import mp


def tanh_sinh_quadrature(f, a, b, eps, max_steps=10, f_derivatives=None):
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

    def fun(t):
        return f(a + 0.5*(b-a) * (t+1))

    fun_derivatives = None
    if f_derivatives:
        fun_derivatives = {
            k: lambda t: fp(a + 0.5*(b-a) * (t+1))
            for k, fp in f_derivatives.items()
            }

    num_digits = int(-mp.log10(eps) + 1)
    mpmath.mp.dps = num_digits

    value_estimates = []
    h = mp.mpf(1)
    success = False
    for level in range(max_steps):
        # For h=1, the error estimate is too optimistic. Hence, start with
        # h=1/2 right away.
        h /= 2

        # We would like to calculate the weights until they are smaller than
        # tau = eps**2, i.e.,
        #
        #     h * pi/2 * cosh(h*j) / cosh(pi/2 * sinh(h*j))**2 < tau.
        #
        # To streamline the computation, j is estimated in advance. The only
        # assumption we're making is that h*j>>1 such that exp(-h*j) can be
        # neglected. With this, the above becomes
        #
        #     tau > h * pi/2 * exp(h*j)/2 / cosh(pi/2 * exp(h*j)/2)
        #
        # and further
        #
        #     tau > h * pi/2 * exp(h*j) / exp(pi/2 * exp(h*j)).
        #
        # Calling z = - pi/2 * exp(h*j), one gets
        #
        #     tau > -h*z / exp(-z)
        #
        # This inequality is fulfilled exactly if z < W(-tau/h) with W being
        # the (-1)-branch of the Lambert-W function IF e*tau < h (which we can
        # assume since `tau` will generally be small). We finally get
        #
        #     j > ln(-W(-tau/h) * 2 / pi) / h.
        #
        assert mp.exp(1)*eps**2 < h
        j = int(mp.ln(-mp.lambertw(-eps**2 / h, -1) * 2 / mp.pi) / h) + 1

        u2 = [mp.pi/2 * mp.sinh(h*jj) for jj in range(j+1)]
        cosh_u2 = [mp.cosh(v) for v in u2]
        weights = [
            h * mp.pi/2 * mp.cosh(h*jj) / v**2
            for jj, v in zip(range(j+1), cosh_u2)
            ]

        # y = 1 - x
        # x = [mp.tanh(v) for v in u2]
        y = [1 / (mp.exp(v) * c) for v, c in zip(u2, cosh_u2)]

        # Perform the integration.
        # The summands are listed such that the points are in ascending order.
        # (The slice expression [-1:0:-1] cuts the first entry and reverses the
        # array.)
        summands = (
            [fun(yy-1) * w for yy, w in zip(y[-1:0:-1], weights[-1:0:-1])]
            + [fun(1-yy) * w for yy, w in zip(y, weights)]
            )
        value_estimates.append(mp.fsum(summands) * (b - a)*0.5)

        # error estimation
        if fun_derivatives:
            error_estimate = _error_estimate1(h, j, f, fun_derivatives)
        else:
            error_estimate = _error_estimate2(
                level, value_estimates, summands, eps
                )

        print(mp.mpf(2)/3 - value_estimates[-1], error_estimate)

        if abs(error_estimate) < eps:
            success = True
            break

    assert success
    return value_estimates[-1], error_estimate


def _error_estimate1(h, j, f, f_derivatives):
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
    return s * h * (h/2/mp.pi)**2


def _error_estimate2(level, value_estimates, summands, eps):
    # "less formal" error estimation after Bailey,
    # <http://www.davidhbailey.com/dhbpapers/dhb-tanh-sinh.pdf>
    if level <= 1:
        error_estimate = 1
    elif value_estimates[0] == value_estimates[-1]:
        error_estimate = 0
    else:
        # d1 = mp.log10(abs(value_estimates[-1] - value_estimates[-2]))
        # d2 = mp.log10(abs(value_estimates[-1] - value_estimates[-3]))
        # d3 = mp.log10(eps * max([abs(x) for x in summands]))
        # d4 = mp.log10(max(abs(summands[0]), abs(summands[-1])))
        # d = max(d1**2 / d2, 2*d1, d3, d4)
        # error_estimate = 10**d
        e1 = abs(value_estimates[-1] - value_estimates[-2])
        e2 = abs(value_estimates[-1] - value_estimates[-3])
        e3 = eps * max([abs(x) for x in summands])
        e4 = max(abs(summands[0]), abs(summands[-1]))
        error_estimate = max(e1**(mp.log(e1)/mp.log(e2)), e1**2, e3, e4)

    return error_estimate
