# -*- coding: utf-8 -*-
#
import math

import numpy
from orthopy import jacobi_recursion_coefficients, scheme_from_coefficients

from .. import helpers
from .gauss_legendre import GaussLegendre


class GaussKronrod(object):
    '''
    Gauss-Kronrod quadrature; see
    <https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula>.

    Besides points and weights, this class provides the weights of the
    corresponding Gauss-Legendre scheme in self.gauss_weights.

    Code adapted from
    <https://www.cs.purdue.edu/archives/2002/wxg/codes/r_kronrod.m>,
    <https://www.cs.purdue.edu/archives/2002/wxg/codes/kronrod.m>.
    See

    Calculation of Gauss-Kronrod quadrature rules,
    Dirk P. Laurie,
    Math. Comp. 66 (1997), 1133-1145,
    <https://doi.org/10.1090/S0025-5718-97-00861-2>

    Abstract:
    The Jacobi matrix of the $(2n+1)$-point Gauss-Kronrod quadrature rule for a
    given measure is calculated efficiently by a five-term recurrence relation.
    The algorithm uses only rational operations and is therefore also useful
    for obtaining the Jacobi-Kronrod matrix analytically. The nodes and weights
    can then be computed directly by standard software for Gaussian quadrature
    formulas.
    '''
    def __init__(self, n, a=0.0, b=0.0):
        # The general scheme is:
        # Get the Jacobi recursion coefficients, get the Kronrod vectors alpha
        # and beta, and hand those off to scheme_from_coefficients. There, the
        # eigenproblem for a tridiagonal matrix with alpha and beta is solved
        # to retrieve the points and weights.
        length = int(math.ceil(3*n/2.0)) + 1
        self.degree = 2*length + 1
        alpha, beta = jacobi_recursion_coefficients(length, a, b)
        a, b = self.r_kronrod(n, alpha, beta)
        x, w = scheme_from_coefficients(a, b)
        # sort by x
        i = numpy.argsort(x)
        self.points = x[i]
        self.weights = w[i]
        return

    # pylint: disable=no-self-use
    def r_kronrod(self, n, a0, b0):
        assert len(a0) == int(math.ceil(3*n/2.0)) + 1
        assert len(b0) == int(math.ceil(3*n/2.0)) + 1

        a = numpy.zeros(2*n+1)
        b = numpy.zeros(2*n+1)

        k = int(math.floor(3*n/2.0)) + 1
        a[:k] = a0[:k]
        k = int(math.ceil(3*n/2.0)) + 1
        b[:k] = b0[:k]
        s = numpy.zeros(int(math.floor(n/2.0)) + 2)
        t = numpy.zeros(int(math.floor(n/2.0)) + 2)
        t[1] = b[n+1]
        for m in range(n-1):
            k0 = int(math.floor((m+1)/2.0))
            k = numpy.arange(k0, -1, -1)
            L = m - k
            s[k+1] = numpy.cumsum(
                (a[k+n+1] - a[L])*t[k+1] + b[k+n+1]*s[k] - b[L]*s[k+1]
                )
            s, t = t, s

        j = int(math.floor(n/2.0)) + 1
        s[1:j+1] = s[:j]
        for m in range(n-1, 2*n-2):
            k0 = m+1-n
            k1 = int(math.floor((m-1)/2.0))
            k = numpy.arange(k0, k1 + 1)
            L = m - k
            j = n-1-L
            s[j+1] = numpy.cumsum(
                -(a[k+n+1] - a[L])*t[j+1] - b[k+n+1]*s[j+1] + b[L]*s[j+2]
                )
            j = j[-1]
            k = int(math.floor((m+1)/2.0))
            if m % 2 == 0:
                a[k+n+1] = a[k] + (s[j+1] - b[k+n+1]*s[j+2]) / t[j+2]
            else:
                b[k+n+1] = s[j+1] / s[j+2]
            s, t = t, s

        a[2*n] = a[n-1] - b[2*n] * s[1]/t[1]
        return a, b


# pylint: disable=too-many-locals
def _gauss_kronrod_integrate(k, f, interval, sumfun=helpers.kahan_sum):
    def _scale_points(points, interval):
        alpha = 0.5 * (interval[1] - interval[0])
        beta = 0.5 * (interval[0] + interval[1])
        return (numpy.multiply.outer(points, alpha) + beta).T

    def _integrate(values, weights, interval_length, sumfun=helpers.kahan_sum):
        '''Integration with point values explicitly specified.
        '''
        out = sumfun(weights * values, axis=-1)
        return 0.5 * interval_length * out

    # Compute the integral estimations according to Gauss and Gauss-Kronrod,
    # sharing the function evaluations
    scheme = GaussKronrod(k)
    gauss_weights = GaussLegendre(k).weights
    point_vals_gk = f(_scale_points(scheme.points, interval))
    point_vals_g = point_vals_gk[..., 1::2]
    alpha = abs(interval[1] - interval[0])
    val_gauss_kronrod = \
        _integrate(point_vals_gk, scheme.weights, alpha, sumfun=sumfun)
    val_gauss = _integrate(point_vals_g, gauss_weights, alpha, sumfun=sumfun)

    # Get an error estimate. According to
    #
    #   A Review of Error Estimation in Adaptive Quadrature
    #   Pedro Gonnet,
    #   ACM Computing Surveys (CSUR) Surveys,
    #   Volume 44, Issue 4, August 2012
    #   <https://doi.org/10.1145/2333112.2333117>,
    #   <https://arxiv.org/pdf/1003.4629.pdf>
    #
    # the classicial QUADPACK still compares favorably with other approaches.
    average = val_gauss_kronrod / alpha
    point_vals_abs = abs(point_vals_gk - average[..., None])
    I_tilde = _integrate(point_vals_abs, scheme.weights, alpha, sumfun=sumfun)
    # The exponent 1.5 is chosen such that (200*x)**1.5 is approximately x at
    # 1.0e-6, the machine precision on IEEE 754 32-bit floating point
    # arithmentic. This could be adapted to
    #
    #   eps = numpy.finfo(float).eps
    #   exponent = numpy.log(eps) / numpy.log(200*eps)
    #
    error_estimate = \
        I_tilde * numpy.minimum(
            numpy.ones(I_tilde.shape),
            (200 * abs(val_gauss_kronrod - val_gauss) / I_tilde)**1.5
            )
    return val_gauss_kronrod, val_gauss, error_estimate
