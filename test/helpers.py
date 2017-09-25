# -*- coding: utf-8 -*-
#
from __future__ import division

import math
import numpy

import scipy.special


def check_degree_1d(
        quadrature, exact, max_degree, tol=1.0e-14
        ):
    val = quadrature(
        lambda x: [x**degree for degree in range(max_degree+1)]
        ).flatten()
    exact_val = numpy.array([exact(degree) for degree in range(max_degree+1)])
    eps = numpy.finfo(float).eps
    # check relative error
    # Allow 1e1 over machine precision.
    alpha = abs(exact_val) * tol + (1e1+tol+exact_val)*eps
    # check where the error is larger than alpha
    is_larger = (exact_val - val) > alpha
    return numpy.where(is_larger)[0] - 1 if any(is_larger) else max_degree


def check_degree(
        quadrature, exact, exponents_creator, max_degree, tol=1.0e-14
        ):
    exponents = numpy.concatenate([
        exponents_creator(degree)
        for degree in range(max_degree+1)
        ])

    exact_vals = numpy.array([exact(k) for k in exponents])

    def fun(x):
        # Evaluate monomials.
        # There's a more complex, faster implementation using matmul, exp, log.
        # However, this only works for strictly positive `x`, and requires some
        # tinkering. See below and
        # <https://stackoverflow.com/a/45421128/353337>.
        return numpy.prod(x[..., None] ** exponents.T[:, None], axis=0).T

    # def fun(x):
    #     # Evaluate many monomials `x^k y^l z^m` at many points. Note that
    #     # `x^k y^l z^m = exp(log(x)*k  + log(y)*l + log(z)*m` , i.e., a
    #     # dot-product, for positive x, y, z. With a correction for points
    #     # with nonpositive components, this is used here.
    #     k = exponents
    #     eps = 1.0e-12
    #     out = \
    #        numpy.exp(numpy.matmul(k, numpy.log(abs(x), where=abs(x) > eps)))
    #     odd_k = numpy.zeros(k.shape, dtype=int)
    #     odd_k[k % 2 == 1] = 1
    #     neg_x = numpy.zeros(x.shape, dtype=int)
    #     neg_x[x < 0.0] = 1
    #     negative_count = numpy.dot(odd_k, neg_x)
    #     out *= (-1)**negative_count
    #     pos_k = numpy.zeros(k.shape, dtype=int)
    #     pos_k[k > 0] = 1
    #     zero_x = numpy.zeros(x.shape, dtype=int)
    #     zero_x[x == 0.0] = 1
    #     zero_count = numpy.dot(pos_k, zero_x)
    #     out[zero_count > 0] = 0.0
    #     return out

    vals = quadrature(fun)

    # check relative error
    # The allowance is quite large here, 1e5 over machine precision.
    # Some tests fail if lowered, though.
    # TODO increase precision
    eps = numpy.finfo(float).eps
    alpha = abs(exact_vals) * tol + (1.0e5+tol+exact_vals)*eps
    is_smaller = abs(exact_vals - vals) < alpha

    if numpy.all(is_smaller):
        return max_degree

    k = numpy.where(numpy.logical_not(is_smaller))[0]
    degree = numpy.sum(exponents[k[0]]) - 1
    return degree


def integrate_monomial_over_standard_simplex(k):
    '''The integrals of monomials over the standard triangle and tetrahedron are
    given by

    \\int_T x_0^k0 * x1^k1 = (k0!*k1!) / (2+k0+k1)!,
    \\int_T x_0^k0 * x1^k1 * x2^k2 = (k0!*k1!*k2!) / (4+k0+k1+k2)!,

    see, e.g.,
    A set of symmetric quadrature rules on triangles and tetrahedra,
    Linbo Zhang, Tao Cui and Hui Liu,
    Journal of Computational Mathematics,
    Vol. 27, No. 1 (January 2009), pp. 89-96,
    <https://www.jstor.org/stable/43693493>.

    See, e.g., <https://math.stackexchange.com/q/207073/36678> for a formula in
    all dimensions.
    '''
    # exp-log to account for large values in numerator and denominator
    return math.exp(
        math.fsum(scipy.special.gammaln(k+1))
        - scipy.special.gammaln(sum(k+1) + 1)
        )


def integrate_monomial_over_enr2(k):
    if numpy.any(k % 2 == 1):
        return 0
    return numpy.prod([math.gamma((kk+1) / 2.0) for kk in k])


def integrate_monomial_over_enr(k):
    if numpy.any(k % 2 == 1):
        return 0
    n = len(k)
    return 2 * math.factorial(sum(k) + n - 1) * numpy.prod([
        math.gamma((kk+1) / 2.0) for kk in k
        ]) / math.gamma((sum(k) + n) / 2)
