# -*- coding: utf-8 -*-
#
import mpmath
from mpmath import mp


# def tanh_sinh_singular(
#         f_left=None,
#         f_right=None,
#         ab=None,
#         eps=None,
#         max_steps=10
#         ):
#     '''Tanh-sinh quadrature for integrands with singularities at the
#     endpoints.
#
#     Background:
#     Since tanh-sinh places its nodes very close to the interval boundaries,
#     it is important to evaluate them with some accuracy. This may become a
#     problem if there are singularities at the boundaries, i.e., if the values
#     change very rapidly. To avoid round-off errors, the user must provide the
#     integrand in terms of distance the left or right hade boundary,
#     respectively. For example, if your integrand is f(t) = sqrt(t / (1-t^2))
#     between 0 and 1, you'll have to provide g(s) = sqrt((1-s) / (2*s - s^2))
#     (where s varies between 0 and 1).
#     If there are singularities on both sides, you can provide both f_left and
#     f_right and it will split the integral in two.
#     '''
#     if f_left and not f_right:
#         value, error = _tanh_sinh(f_left, eps, max_steps=max_steps)
#     elif f_right and not f_left:
#         value, error = _tanh_sinh(f_right, eps, max_steps=max_steps)
#     else:
#         assert f_left and f_right
#         value0, error0 = \
#             _tanh_sinh(f_left, 0, ab/2, eps, max_steps=max_steps)
#         value1, error1 = \
#             _tanh_sinh(f_right, 0, ab/2, eps, max_steps=max_steps)
#         value = value0 + value1
#         error = error0 + error1
#     return value, error


def tanh_sinh_regular(f, a, b, eps, max_steps=10, f_derivatives=None):
    if f_derivatives is None:
        f_derivatives = {}

    ba2 = (b-a) / mp.mpf(2)

    f_left = {
        0: lambda s: ba2**0 * f(a + s*ba2),
        1: lambda s: ba2**1 * f_derivatives[1](a + s*ba2),
        2: lambda s: ba2**2 * f_derivatives[2](a + s*ba2),
        }
    f_right = {
        0: lambda s: (-ba2)**0 * f(b - s*ba2),
        1: lambda s: (-ba2)**1 * f_derivatives[1](b - s*ba2),
        2: lambda s: (-ba2)**2 * f_derivatives[2](b - s*ba2),
        }

    value_estimate, error_estimate = _tanh_sinh(
        f_left, f_right, eps,
        max_steps=max_steps
        )
    return value_estimate * ba2, error_estimate * ba2


# pylint: disable=too-many-arguments, too-many-locals
def _tanh_sinh(f_left, f_right, eps, max_steps=10):
    '''Integrate a function `f` between `a` and `b` with accuracy `eps`. The
    function `f` is given in terms of two functions `f_left` and `f_right`
    where `f_left(s) = f(a - s*(a-b)/2)`, i.e., `f` linearly scaled such that
    `f_left(0) = a`, `f_left(2) = b` (`f_right` likewise).

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
        #     tau > h * pi/2 * exp(h*j)/2 / cosh(pi/2 * exp(h*j)/2)**2
        #
        # and further
        #
        #     tau > h * pi * exp(h*j) / exp(pi/2 * exp(h*j)).
        #
        # Calling z = - pi/2 * exp(h*j), one gets
        #
        #     tau > -2*h*z * exp(z)
        #
        # This inequality is fulfilled exactly if z = W(-tau/h/2) with W being
        # the (-1)-branch of the Lambert-W function IF exp(1)*tau < 2*h (which
        # we can assume since `tau` will generally be small). We finally get
        #
        #     j > ln(-2/pi * W(-tau/h/2)) / h.
        #
        assert mp.exp(1)*eps**2 < 2*h
        j = int(mp.ln(-2/mp.pi * mp.lambertw(-eps**2/h/2, -1)) / h) + 1

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
            [f_left[0](yy) * w for yy, w in zip(y[-1:0:-1], weights[-1:0:-1])]
            + [f_right[0](yy) * w for yy, w in zip(y, weights)]
            )
        value_estimates.append(mp.fsum(summands))

        # error estimation
        if 1 in f_left and 2 in f_left:
            assert 1 in f_right and 2 in f_right
            error_estimate = _error_estimate1(h, j, f_left, f_right)
        else:
            error_estimate = _error_estimate2(
                level, value_estimates, summands, eps
                )

        # exact = (mp.exp(mp.pi/2) - 1)/2
        exact = mp.pi/4
        print(value_estimates[-1] * 0.5, exact)
        print(value_estimates[-1] * 0.5 - exact, error_estimate * 0.5)
        print
        if level == 2:
            exit(1)

        if abs(error_estimate) < eps:
            success = True
            break

    assert success
    return value_estimates[-1], error_estimate


def _error_estimate1(h, j, f_left, f_right):
    # Pretty accurate error estimation:
    #
    #   E(h) = h * (h/2/pi)**2 * sum_{-N}^{+N} F''(h*j)
    #
    # with
    #
    #   F(t) = f(g(t)) * g'(t),
    #   g(t) = tanh(pi/2 sinh(t)).
    #

    # def g(t):
    #     return mp.tanh(mp.pi/2 * mp.sinh(t))
    # y = 1 - g(t)
    def y(t):
        u2 = mp.pi/2 * mp.sinh(t)
        return 1 / (mp.exp(u2) * mp.cosh(u2))

    def dy_dt(t):
        return -mp.pi/2 * mp.cosh(t) / mp.cosh(mp.pi/2 * mp.sinh(t))**2

    def d2y_dt2(t):
        return -mp.pi/2 * (
            + mp.sinh(t)
            - mp.pi * mp.cosh(t)**2 * mp.tanh(mp.pi/2 * mp.sinh(t))
            ) / mp.cosh(mp.pi/2 * mp.sinh(t))**2

    def d3y_dt3(t):
        sinh_sinh = mp.sinh(mp.pi/2 * mp.sinh(t))
        cosh_sinh = mp.cosh(mp.pi/2 * mp.sinh(t))
        tanh_sinh = mp.tanh(mp.pi/2 * mp.sinh(t))
        return -mp.pi/4 * mp.cosh(t) * (
            + 2 * cosh_sinh
            - 2 * mp.pi**2 * mp.cosh(t)**2 / cosh_sinh
            + mp.pi**2 * mp.cosh(t)**2 * cosh_sinh
            + mp.pi**2 * mp.cosh(t)**2 * tanh_sinh * sinh_sinh
            - 6 * mp.pi * mp.sinh(t) * sinh_sinh
            ) / cosh_sinh**3

    # TODO reuse yt, g*
    def F2(f, t):
        '''Second derivative of F(t) = f(g(t)) * g'(t).
        '''
        yt = y(t)
        y1 = dy_dt(t)
        y2 = d2y_dt2(t)
        y3 = d3y_dt3(t)
        return y1**3 * f[2](yt) + 3*y1*y2 * f[1](yt) + y3 * f[0](yt)

    t = [h * jj for jj in range(j+1)]
    summands = [F2(f_left, tt) for tt in t[1:]] + [F2(f_right, tt) for tt in t]
    return h * (h/2/mp.pi)**2 * mp.fsum(summands)


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
