# -*- coding: utf-8 -*-
#
import math
import numpy


def _jacobi_recursion_coefficients(n, a, b):
    '''
    Generate the recursion coefficients alpha_k, beta_k

    P_{k+1}(x) = (x-alpha_k)*P_{k}(x) - beta_k P_{k-1}(x)

    for the Jacobi polynomials which are orthogonal on [-1, 1]
    with respect to the weight w(x)=[(1-x)^a]*[(1+x)^b].

    Adapted from the MATLAB code by Dirk Laurie and Walter Gautschi
    http://www.cs.purdue.edu/archives/2002/wxg/codes/r_jacobi.m
    and from Greg van Winckel's
    https://github.com/gregvw/orthopoly-quadrature/blob/master/rec_jacobi.pyx
    '''
    assert a > -1.0 or b > -1.0
    assert n >= 1

    mu = 2.0**(a+b+1.0) \
        * numpy.exp(
            math.lgamma(a+1.0) + math.lgamma(b+1.0) - math.lgamma(a+b+2.0)
            )
    nu = (b-a) / (a+b+2.0)

    if n == 1:
        return nu, mu

    N = numpy.arange(1, n)

    nab = 2.0*N + a + b
    alpha = numpy.hstack([nu, (b**2 - a**2) / (nab * (nab + 2.0))])
    N = N[1:]
    nab = nab[1:]
    B1 = 4.0 * (a+1.0) * (b+1.0) / ((a+b+2.0)**2.0 * (a+b+3.0))
    B = 4.0 * (N+a) * (N+b) * N * (N+a+b) / (nab**2.0 * (nab+1.0) * (nab-1.0))
    beta = numpy.hstack((mu, B1, B))
    return alpha, beta


def _gauss(alpha, beta):
    '''
    Compute the Gauss nodes and weights from the recursion coefficients
    associated with a set of orthogonal polynomials.

    Algorithm 726: ORTHPOL–a package of routines for generating orthogonal
    polynomials and Gauss-type quadrature rules,
    W. Gautschi,
    ACM Transactions on Mathematical Software (TOMS),
    Volume 20, Issue 1, March 1994,
    Pages 21-62,
    <http://doi.org/10.1145/174603.174605>,
    <http://www.cs.purdue.edu/archives/2002/wxg/codes/gauss.m>,

    and

    http://www.scientificpython.net/pyblog/radau-quadrature
    '''
    from scipy.linalg import eig_banded
    import scipy
    A = numpy.vstack((numpy.sqrt(beta), alpha))
    x, V = eig_banded(A, lower=False)
    w = beta[0]*scipy.real(scipy.power(V[0, :], 2))
    return x, w
