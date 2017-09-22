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

    # What's a good initial h?
    # The larger `h` is chosen, the fewer points will be part of the
    # evaluation. However, we don't want to choose the step size too large
    # since that means less accuracy for the quadrature overall. The idea would
    # then be too choose `h` such that it is just large enough for the first
    # tanh-sinh-step to contain only one point, the midpoint. The expression
    #
    #    j = mp.ln(-2/mp.pi * mp.lambertw(-tau/h/2, -1)) / h
    #
    # hence needs to just smaller than 1. One gets
    #
    #    0 = pi/2 * exp(h) - h - ln(h) - ln(pi/tau)
    #
    # for which there is no analytic solution of this equation, but one can
    # approximate it. Since pi/2 * exp(h) >> h >> ln(h) (for `h` large enough),
    # one can either forget about both h and ln(h) to get
    #
    #     h0 = ln(2/pi * ln(pi/tau))
    #
    # or just scratch ln(h) to get
    #
    #     h1 = ln(tau/pi) - W_{-1}(-tau/2).
    #
    # Both of these suggestions underestimate and `j` will be too large. An
    # approximation that overestimates is obtained by replacing `ln(h)` by `h`,
    #
    #     h2 = 1/2 - log(sqrt(pi/tau)) - W_{-1}(-sqrt(exp(1)*pi*tau) / 4).
    #
    # Application of Newton's method will improve all of these approximations
    # and will also always overestimate such that `j` won't exceed 1 in the
    # first step. Nice!
    h = _solve_expx_x_logx(eps**2, tol=1.0e-10)

    success = False
    for level in range(max_steps+1):
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
        # We do require j to be positive, so -2/pi * W(-tau/h/2) > 1. This
        # translates to the slightly stricter requirement
        #
        #     tau * exp(pi/2) < pi * h,
        #
        # i.e., h needs to be about 1.531 times larger than tau (not only 1.359
        # times as the previous bound suggested).
        #
        # Note further that h*j is ever decreasing as h decreases.
        assert eps**2 * mp.exp(mp.pi/2) < mp.pi*h
        j = int(mp.ln(-2/mp.pi * mp.lambertw(-eps**2/h/2, -1)) / h)

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

        if level == 0:
            # The root level only contains one node, the midpoint; function
            # values of f_left and f_right are equal here.
            assert len(weights) == 1
            assert len(y0) == 1
            value_estimates = [weights[0] * f_left[0](y0[0])]
            error_estimate = 1
        else:
            # Perform the integration.
            # The summands are listed such that the points are in ascending
            # order. (The slice expression [-1:0:-1] cuts the first entry and
            # reverses the array.)
            fly = numpy.array([f_left[0](yy) for yy in y0])
            fry = numpy.array([f_right[0](yy) for yy in y0])
            summands = numpy.concatenate([
                fly[1::2] * weights[1::2], fry[1::2] * weights[1::2]
                ])
            value_estimates.append(
                value_estimates[-1]/2 + mp.fsum(summands)
                )

            # error estimation
            if 1 in f_left and 2 in f_left:
                assert 1 in f_right and 2 in f_right
                last_error_estimate = error_estimate
                error_estimate = _error_estimate1(
                    h, sinh_t, cosh_t, cosh_sinh_t, y0, y1,
                    fly, fry, f_left, f_right, alpha, last_error_estimate
                    )
            else:
                index_leftmost = len(weights) - 2
                index_rightmost = -1
                error_estimate = _error_estimate2(
                    eps, value_estimates,
                    summands, index_leftmost, index_rightmost
                    )

        if abs(error_estimate) < eps:
            success = True
            break

        h /= 2

    assert success
    return value_estimates[-1], error_estimate


def _error_estimate1(
        h, sinh_t, cosh_t, cosh_sinh_t, y0, y1,
        fly, fry, f_left, f_right, alpha, last_estimate
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

    fl1_y = numpy.array([f_left[1](yy) for yy in y0])
    fl2_y = numpy.array([f_left[2](yy) for yy in y0])

    fr1_y = numpy.array([f_right[1](yy) for yy in y0])
    fr2_y = numpy.array([f_right[2](yy) for yy in y0])

    # Second derivative of F(t) = f(g(t)) * g'(t).
    summands = numpy.concatenate([
        (y3 * fly + 3*y1*y2 * fl1_y + y1**3 * fl2_y)[1:],
        y3 * fry + 3*y1*y2 * fr1_y + y1**3 * fr2_y,
        ])
    return h * (h/2/mp.pi)**2 * mp.fsum(summands)


def _error_estimate2(
        eps, value_estimates,
        summands, index_leftmost, index_rightmost
        ):
    # "less formal" error estimation after Bailey,
    # <http://www.davidhbailey.com/dhbpapers/dhb-tanh-sinh.pdf>
    if len(value_estimates) < 3:
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
        e4 = max(abs(summands[index_leftmost]), abs(summands[index_rightmost]))
        error_estimate = max(e1**(mp.log(e1)/mp.log(e2)), e1**2, e3, e4)

    return error_estimate


def _solve_expx_x_logx(tau, tol, max_steps=10):
    '''Solves the equation

    log(pi/tau) = pi/2 * exp(x) - x - log(x)

    approximately using Newton's method. The approximate solution is guaranteed
    to overestimate.
    '''
    x = mp.log(2/mp.pi * mp.log(mp.pi/tau))
    # x = mp.log(tau/mp.pi) -  mp.lambertw(-tau/2, -1))
    # x = mp.mpf(1)/2 \
    #    - mp.log(mp.sqrt(mp.pi/tau)) \
    #    - mp.lambertw(-mp.sqrt(mp.exp(1)*mp.pi*tau)/4, -1)

    def f0(x):
        return mp.pi/2 * mp.exp(x) - x - mp.log(x*mp.pi/tau)

    def f1(x):
        return mp.pi/2 * mp.exp(x) - 1 - mp.mpf(1)/x

    f0x = f0(x)
    success = False
    # At least one step is performed. This is required for the guarantee of
    # overestimation.
    for k in range(max_steps):
        x -= f0x / f1(x)
        f0x = f0(x)
        if abs(f0x) < tol:
            success = True
            break

    assert success
    return x
