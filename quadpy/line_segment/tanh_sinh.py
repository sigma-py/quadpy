# -*- coding: utf-8 -*-
#
from mpmath import mp
import numpy


# pylint: disable=too-many-arguments
def tanh_sinh_quadrature(f, a, b, eps, max_steps=10, f_derivatives=None):
    if f_derivatives is None:
        f_derivatives = {}

    f_left = {
        0: lambda s: f(a + s),
        }
    if 1 in f_derivatives:
        f_left[1] = lambda s: f_derivatives[1](a + s)
    if 2 in f_derivatives:
        f_left[2] = lambda s: f_derivatives[2](a + s)

    f_right = {
        0: lambda s: f(b - s),
        }
    if 1 in f_derivatives:
        f_right[1] = lambda s: -f_derivatives[1](b - s)
    if 2 in f_derivatives:
        f_right[2] = lambda s: +f_derivatives[2](b - s)

    value_estimate, error_estimate = tanh_sinh_lr(
        f_left, f_right, b-a, eps,
        max_steps=max_steps
        )
    return value_estimate, error_estimate


# pylint: disable=too-many-arguments, too-many-locals
def tanh_sinh_lr(f_left, f_right, alpha, eps, max_steps=10):
    '''Integrate a function `f` between `a` and `b` with accuracy `eps`. The
    function `f` is given in terms of two functions

        * `f_left(s) = f(a + s)`, i.e., `f` linearly scaled such that
          `f_left(0) = a`, `f_left(b-a) = b`,

        * `f_right(s) = f(b - s)`, i.e., `f` linearly scaled such that
          `f_right(0) = b`, `f_left(b-a) = a`.

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
    mp.dps = num_digits

    alpha2 = alpha / mp.mpf(2)

    value_estimates = []
    h = mp.mpf(1)
    success = False
    for level in range(max_steps):
        print(level)
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

        t = numpy.array([h * jj for jj in range(j+1)])

        sinh_t = mp.pi/2 * numpy.array(map(mp.sinh, t))
        cosh_t = mp.pi/2 * numpy.array(map(mp.cosh, t))
        cosh_sinh_t = numpy.array(map(mp.cosh, sinh_t))

        # y = alpha/2 * (1 - x)
        # x = [mp.tanh(v) for v in u2]
        exp_sinh_t = numpy.array(map(mp.exp, sinh_t))

        y0 = alpha2 / exp_sinh_t / cosh_sinh_t
        y1 = -alpha2 * cosh_t / cosh_sinh_t**2

        weights = -h * y1

        # Perform the integration.
        # The summands are listed such that the points are in ascending order.
        # (The slice expression [-1:0:-1] cuts the first entry and reverses the
        # array.)
        fly = numpy.array([f_left[0](yy) for yy in y0])
        fry = numpy.array([f_right[0](yy) for yy in y0])
        summands = numpy.concatenate([
            fly[-1:0:-1] * weights[-1:0:-1],
            fry*weights
            ])
        value_estimates.append(mp.fsum(summands))

        # error estimation
        if 1 in f_left and 2 in f_left:
            assert 1 in f_right and 2 in f_right
            error_estimate = _error_estimate1(
                h, t, sinh_t, cosh_t, cosh_sinh_t, y0, y1,
                fly, fry, f_left, f_right, alpha
                )
        else:
            error_estimate = _error_estimate2(
                level, value_estimates, summands, eps
                )

        if abs(error_estimate) < eps:
            success = True
            break

    assert success
    return value_estimates[-1], error_estimate


def _error_estimate1(
        h, t, sinh_t, cosh_t, cosh_sinh_t, y0, y1,
        fly, fry, f_left, f_right, alpha
        ):
    '''
    A pretty accurate error estimation is

      E(h) = h * (h/2/pi)**2 * sum_{-N}^{+N} F''(h*j)

    with

      F(t) = f(g(t)) * g'(t),
      g(t) = tanh(pi/2 sinh(t)).
    '''
    alpha2 = alpha / mp.mpf(2)

    sinh_sinh_t = numpy.array(map(mp.sinh, sinh_t))
    tanh_sinh_t = sinh_sinh_t / cosh_sinh_t

    # More derivatives of y = 1-g(t).
    y2 = -alpha2 * (sinh_t - 2 * cosh_t**2 * tanh_sinh_t) / cosh_sinh_t**2
    y3 = -alpha2 * cosh_t * (
        + cosh_sinh_t
        - 4 * cosh_t**2 / cosh_sinh_t
        + 2 * cosh_t**2 * cosh_sinh_t
        + 2 * cosh_t**2 * tanh_sinh_t * sinh_sinh_t
        - 6 * sinh_t * sinh_sinh_t
        ) / cosh_sinh_t**3

    def F2(f, f0_y, y):
        '''Second derivative of F(t) = f(g(t)) * g'(t).
        '''
        return (
            + y[3] * f0_y
            + 3*y[1]*y[2] * f[1](y[0])
            + y[1]**3 * f[2](y[0])
            )

    y = numpy.array([y0, y1, y2, y3]).T

    # fl1_y = numpy.array([f_left[1](yy) for yy in y])
    # fl2_y = numpy.array([f_left[2](yy) for yy in y])
    # #
    # fr1_y = numpy.array([f_right[1](yy) for yy in y])
    # fr2_y = numpy.array([f_right[2](yy) for yy in y])
    summands = (
        [F2(f_left, f0, yy) for f0, yy in zip(fly[1:], y[1:])]
        + [F2(f_right, f0, yy) for f0, yy in zip(fry, y)]
        )
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
