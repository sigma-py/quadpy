# -*- coding: utf-8 -*-
#
from __future__ import division

import math
import numpy


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


def evaluate_all_monomials(x, max_degree):
    '''Evaluate all polynomials of up to degree `max_degree` at point(s) `x`.
    '''
    def rec(exponents, vals, x):
        if len(exponents) == 0 or len(exponents[0]) == 0:
            return [], []

        idx_leading_zero = [
            k for k in range(len(exponents)) if exponents[k][0] == 0
            ]
        exponents_with_leading_zero = \
            [exponents[k][1:] for k in idx_leading_zero]
        val1 = vals[idx_leading_zero]
        x1 = x[1:]
        out1, vals1 = rec(exponents_with_leading_zero, val1, x1)

        # increment leading exponent by 1
        out = [[e[0]+1] + e[1:] for e in exponents]
        vals0 = vals * x[0]

        out += [[0] + e for e in out1]
        out_vals = numpy.concatenate([vals0, vals1])
        return out, out_vals

    dim = x.shape[0]

    # Initialization, level 0
    exponents = [dim * [0]]
    vals = numpy.array(numpy.ones(x.shape[1:]))

    all_vals = []
    all_exponents = []

    all_exponents.append(exponents)
    all_vals.append(vals)
    for k in range(max_degree):
        exponents, vals = rec(exponents, vals, x)
        all_exponents.append(exponents)
        all_vals.append(vals)

    return all_vals


def check_degree(
        quadrature, exact, exponents_creator, dim, max_degree, tol=1.0e-14
        ):
    exponents = numpy.concatenate([
        exponents_creator(degree)
        for degree in range(max_degree+1)
        ])
    exact_vals = numpy.array([exact(k) for k in exponents])

    # def fun(x):
    #   # Evaluate monomials.
    #   # There's a more complex, faster implementation using matmul, exp, log.
    #   # However, this only works for strictly positive `x`, and requires some
    #   # tinkering. See below and
    #   # <https://stackoverflow.com/a/45421128/353337>.
    #   return numpy.prod(x[..., None] ** exponents.T[:, None], axis=0).T

    vals = quadrature(evaluate_all_monomials)

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
