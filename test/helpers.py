# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy


def check_degree_1d(
        quadrature, exact, exponents_creator, max_degree, tol=1.0e-14
        ):
    for degree in range(max_degree+1):
        val = quadrature(lambda x: x**degree)
        exact_val = exact(degree)
        # check relative error
        eps = numpy.finfo(float).eps
        # Allow 1e2 over machine precision.
        alpha = abs(exact_val) * tol + (1e2+tol+exact_val)*eps
        if abs(exact_val - val) > alpha:
            return degree - 1
    return max_degree


def check_degree(
        quadrature, exact, exponents_creator, max_degree, tol=1.0e-14
        ):
    for degree in range(max_degree+1):
        for k in exponents_creator(degree):
            val = quadrature(
                lambda x: sympy.prod([x[i]**k[i] for i in range(len(k))])
                )
            exact_val = exact(k)
            # check relative error
            # The allowance is quite large here, 1e5 over machine precision.
            # Some test fail if lowered, though.
            # TODO increase precision
            eps = numpy.finfo(float).eps
            alpha = abs(exact_val) * tol + (1.0e5+tol+exact_val)*eps
            if abs(exact_val - val) > alpha:
                return degree - 1
    return max_degree


def create_monomial_exponents2(degree):
    '''Returns a list of all monomials exponents of degree :degree:.
    '''
    return [(degree-k, k) for k in range(degree+1)]


def create_monomial_exponents3(degree):
    '''Returns a list of all monomial exponents of degree :degree:.
    '''
    return numpy.array([
        [degree-i-j, i, j]
        for i in range(degree+1)
        for j in range(degree-i+1)
        ],
        dtype=int
        )


def integrate_monomial_over_unit_circle(k):
    '''The integral

    I = \int_0^2pi cos(phi)**k[0] sin(phi)**k[1]

    equals 0 if any of k[0], k[1] is odd. If both are even, we make use of

    I = 4 \int_0^pi/2 cos(phi)**k[0] sin(phi)**k[1]
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
