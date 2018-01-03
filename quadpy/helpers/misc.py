# -*- coding: utf-8 -*-
#
from __future__ import division

import math
import numpy
import scipy.special
import sympy


def untangle(data):
    weights, points = zip(*data)
    return (
        numpy.concatenate(points),
        numpy.repeat(weights, [len(grp) for grp in points])
        )


def n_outer(a):
    '''Given a list (tuple, array) of arrays, this method computes their outer
    product. If the dimension of the input arrays is larger than one, the
    product is formed across the first dimension; all other dimensions must
    coincide in size.

    Examples:
    n_outer([np.ones(4), np.ones(5)]).shape == (4, 5)
    n_outer([np.ones(4), np.ones(5), np.ones(6)]).shape == (4, 5, 6)
    n_outer([np.ones(4, 3, 7), np.ones(5, 3, 7)]).shape == (4, 5, 3, 7)
    '''
    # <https://stackoverflow.com/a/45376730/353337>
    d = len(a)

    # If the elements are more than one-dimensional, assert that the extra
    # dimensions are all equal.
    s0 = a[0].shape
    for arr in a:
        assert s0[1:] == arr.shape[1:]

    out = a[0]
    for k in range(1, d):
        # Basically outer products. Checkout `numpy.outer`'s implementation for
        # comparison.
        out = numpy.multiply(
                # Insert a newaxis after k `:`
                out[(slice(None),) * k + (numpy.newaxis,)],
                # Insert a newaxis at the beginning
                a[k][numpy.newaxis],
                )
    return out


def kahan_sum(a, axis=0):
    '''Kahan summation of the numpy array `a` along axis `axis`.
    '''
    # See <https://en.wikipedia.org/wiki/Kahan_summation_algorithm> for
    # details.
    k = axis % len(a.shape)
    s = numpy.zeros(a.shape[:axis] + a.shape[k+1:])
    c = numpy.zeros(s.shape)
    for i in range(a.shape[axis]):
        # http://stackoverflow.com/a/42817610/353337
        y = a[(slice(None),) * k + (i,)] - c
        t = s + y
        c = (t - s) - y
        s = t.copy()
    return s


def kahan_dot(a, b):
    '''The dot product performed as a Kahan sum.
    '''
    return kahan_sum(numpy.moveaxis(b * a, -1, 0))


# pylint: disable=too-many-arguments, too-many-locals
def compute_dobrodeev(n, I0, I2, I22, I4, pm_type, i, j, k, symbolic=False):
    '''Compute some helper quantities used in

    L.N. Dobrodeev,
    Cubature rules with equal coefficients for integrating functions with
    respect to symmetric domains,
    USSR Computational Mathematics and Mathematical Physics,
    Volume 18, Issue 4, 1978, Pages 27-34,
    <https://doi.org/10.1016/0041-5553(78)90064-2>.
    '''
    t = 1 if pm_type == 'I' else -1

    binomial = sympy.binomial if symbolic else scipy.special.binom
    fact = sympy.factorial if symbolic else math.factorial
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    L = binomial(n, i) * 2**i
    M = fact(n) // (fact(j) * fact(k) * fact(n-j-k)) * 2**(j+k)
    N = L + M
    F = I22/I0 - I2**2/I0**2 + (I4/I0 - I22/I0) / n
    R = -(j+k-i) / i * I2**2/I0**2 + (j+k-1)/n * I4/I0 - (n-1)/n * I22/I0
    H = 1/i * (
        (j+k-i) * I2**2/I0**2 + (j+k)/n * ((i-1) * I4/I0 - (n-1)*I22/I0)
        )
    Q = L/M*R + H - t * 2*I2/I0 * (j+k-i)/i * sqrt(L/M*F)

    G = 1/N
    a = sqrt(n/i * (I2/I0 + t * sqrt(M/L*F)))
    b = sqrt(n/(j+k) * (I2/I0 - t * sqrt(L/M*F) + t * sqrt(k/j*Q)))
    c = sqrt(n/(j+k) * (I2/I0 - t * sqrt(L/M*F) - t * sqrt(j/k*Q)))
    return G, a, b, c
