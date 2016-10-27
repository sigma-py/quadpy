# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy


def check_degree(quadrature, exact, max_degree, tol=1.0e-10):
    for degree in range(max_degree+1):
        for k in create_monomial_exponents2(degree):
            val = quadrature(lambda x: x[0]**k[0] * x[1]**k[1])
            exact_val = exact(k)
            if abs(exact_val - val) > 1.0e-10:
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

    equals 0 if any of k[0], k[1] is odd. If both are even, then

    I = \int_0^2pi sin**2m * cos**2n
      = \int_0^2pi (1-cos**2)**m * cos**2n
      = \sum_k=0^m (m over k) \int cos**2(k+n).

    With
      \int cos**n = (n-1)/n \int cos**(n-2)
    one gets
        \int_0^2pi cos**n
      = (n-1)!! / n!! * 2*pi
      = (2n)! / (2**n * n!)**2 * 2*pi,
    so
      I = \sum_k=0^m (-1)^i (m over k) (2(k+n))!
                     / (2**(k+n) * (k+n)!)**2 * 2*pi

    The quotient can be computed via log-gamma, i.e.,

      (2d)! / (2**d * d!)**2
    = exp(lgamma(2d) - 2d log(2) - 2 lgamma(d)).
    '''
    if any(numpy.array(k) % 2 == 1):
        return 0.0

    m = k[0] // 2
    n = k[1] // 2
    return 2*math.pi * sum([
            (-1)**i
            * sympy.binomial(m, i)
            * math.factorial(2*(i+n))
            / (2**(i+n) * math.factorial(i+n))**2
            for i in range(m+1)
            ])
    # return 2*math.pi * sum([
    #         (-1)**i
    #         * sympy.binomial(m, i)
    #         * math.exp(
    #             + math.lgamma(2*(i+n))
    #             - 2*math.lgamma(i+n)
    #             - 2*(i+n)*math.log(2.0)
    #             )
    #         ] for i in range(m+1)
    #         )
