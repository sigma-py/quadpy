# -*- coding: utf-8 -*-
#
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


def check_degree(
        quadrature, exact, exponents_creator, max_degree, tol=1.0e-14
        ):
    exponents = numpy.concatenate([
        exponents_creator(degree)
        for degree in range(max_degree+1)
        ])

    exact_vals = numpy.array([exact(k) for k in exponents])

    # def fun(x):
    #     # Naive evaluation of the monomials.
    #     return numpy.prod([
    #         numpy.power.outer(x[i], exponents[:, i]) for i in range(len(x))
    #         ], axis=0).T

    def fun(x):
        # <https://stackoverflow.com/a/45421128/353337>
        # Evaluate many monomials `x^k y^l z^m` at many points. Note that `x^k
        # y^l z^m = exp(log(x)*k  + log(y)*l + log(z)*m` , i.e., a dot-product,
        # for positive x, y, z. With a correction for points with nonpositive
        # components, this is used here.
        k = exponents

        out = numpy.exp(numpy.matmul(k, numpy.log(abs(x), where=abs(x) > 0.0)))

        odd_k = numpy.zeros(k.shape, dtype=int)
        odd_k[k % 2 == 1] = 1
        neg_x = numpy.zeros(x.shape, dtype=int)
        neg_x[x < 0.0] = 1
        negative_count = numpy.dot(odd_k, neg_x)
        out *= (-1)**negative_count

        pos_k = numpy.zeros(k.shape, dtype=int)
        pos_k[k > 0] = 1
        zero_x = numpy.zeros(x.shape, dtype=int)
        zero_x[x == 0.0] = 1
        zero_count = numpy.dot(pos_k, zero_x)
        out[zero_count > 0] = 0.0
        return out

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


def partition(n, d):
    '''Create all nonnegative tuples of length d which sum up to n.
    '''
    # <https://stackoverflow.com/a/45348441/353337>
    def rec(n, d, depth=0):
        if d == depth:
            return [[]]
        return [
            item + [i]
            for i in range(n+1)
            for item in rec(n-i, d, depth=depth+1)
            ]
    return [[n-sum(p)] + p for p in rec(n, d-1)]


def integrate_monomial_over_unit_circle(k):
    '''The integral

    I = \\int_0^2pi \\cos(phi)**k[0] \\sin(phi)**k[1]

    equals 0 if any of k[0], k[1] is odd. If both are even, we make use of

    I = 4 \\int_0^pi/2 \\cos(phi)**k[0] \\sin(phi)**k[1]
      = 2 B(0.5*(k[0]+1), 0.5*(k[1]+1))

    with B(x, y) being the Beta function. It has the representation

    B(x, y) = Gamma(x) Gamma(y) / Gamma(x + y).
    '''
    if any(numpy.array(k) % 2 == 1):
        return 0.0

    return 2 * math.exp(
        + math.lgamma(0.5 * (k[0] + 1))
        + math.lgamma(0.5 * (k[1] + 1))
        - math.lgamma(0.5 * (k[0] + k[1]) + 1)
        )
