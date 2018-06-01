# -*- coding: utf-8 -*-
#
from __future__ import division

from math import factorial as fact

import numpy
import sympy

from .helpers import integrate_monomial_over_unit_nsphere
from ..helpers import untangle, fsd


class Dobrodeev1978(object):
    '''
    L.N. Dobrodeev,
    Cubature rules with equal coefficients for integrating functions with
    respect to symmetric domains,
    USSR Computational Mathematics and Mathematical Physics,
    Volume 18, Issue 4, 1978, Pages 27-34,
    <https://doi.org/10.1016/0041-5553(78)90064-2>.

    Abstract:
    Three-parameter cubature sums, whose values have to be the same as the
    values of the computed integral on all polynomials of degree p <= 5, are
    investigated. The parameters are selected on the basis of a specific
    optimality criterion for a number of symmetric domains of integration, of
    dimensionalities 2 <= n <= 20.
    '''

    def __init__(self, n, symbolic=False):
        # from sympy import sqrt, factorial as fact, Rational as fr

        assert 2 <= n <= 20

        self.name = 'Dobrodeev1978'
        self.degree = 5
        self.dim = n

        dim_config = {
            2: ('I', None, 1, 1),
            3: ('I', None, 1, 1),
            4: (None, 2, None, None),
            5: ('I', 5, 1, 1),
            6: ('I', 5, 1, 1),
            7: (None, 3, None, None),
            8: ('I', 8, 1, 1),
            9: ('I', 9, 1, 1),
            10: ('I', 10, 1, 1),
            11: ('I', 11, 1, 1),
            12: ('I', 5, 1, 1),
            13: ('I', 13, 1, 2),
            14: ('I', 14, 1, 2),
            15: ('I', 15, 1, 2),
            16: ('I', 16, 1, 2),
            17: ('I', 17, 1, 2),
            18: ('I', 18, 1, 3),
            19: ('I', 19, 1, 3),
            20: ('I', 20, 1, 3),
            }

        pm_type, i, j, k = dim_config[n]
        I0 = integrate_monomial_over_unit_nsphere(n * [0], symbolic=symbolic)
        if i is None:
            G, b, c = _generate_jk(n, pm_type, j, k, symbolic)
            data = [(G, fsd(n, (b, j), (c, k)))]
        elif j is None:
            assert k is None
            assert pm_type is None
            G, a = _generate_i(n, i, symbolic)
            data = [(G, fsd(n, (a, i)))]
        else:
            I2 = integrate_monomial_over_unit_nsphere(
                    [2] + (n-1) * [0],
                    symbolic=symbolic
                    )
            I22 = integrate_monomial_over_unit_nsphere(
                    [2, 2] + (n-2) * [0],
                    symbolic=symbolic
                    )
            I4 = integrate_monomial_over_unit_nsphere(
                    [4] + (n-1) * [0],
                    symbolic=symbolic
                    )

            G, a, b, c = _compute_dobrodeev(
                    n, I0, I2, I22, I4, pm_type, i, j, k, symbolic=symbolic
                    )

            data = [
                (G, fsd(n, (a, i))),
                (G, fsd(n, (b, j), (c, k))),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= I0
        return


def _generate_i(n, i, symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    L = fact(n) // fact(i) // fact(n-i) * 2**i
    G = frac(1, L)
    a = sqrt(frac(3, n+2))
    return G, a


def _generate_jk(n, pm_type, j, k, symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    M = fact(n) // fact(j) // fact(k) // fact(n-j-k) * 2**(j+k)
    G = frac(1, M)

    t = 1 if pm_type == 'I' else -1
    b = sqrt(frac(1, j+k) * (1 + t * (k/j * sqrt(3*(j+k)/(n+2) - 1))))
    c = sqrt(frac(1, j+k) * (1 - t * (j/k * sqrt(3*(j+k)/(n+2) - 1))))
    return G, b, c



def _compute_dobrodeev(n, I0, I2, I22, I4, pm_type, i, j, k, symbolic):
    '''Same as the helper function in ..helpers, making use of the fact that
    `F == 0` for the sphere
    '''
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    # TODO prove F==0 analytically
    t = 1 if pm_type == 'I' else -1

    L = fact(n) // (fact(i) * fact(n-i)) * 2**i
    M = fact(n) // (fact(j) * fact(k) * fact(n-j-k)) * 2**(j+k)
    N = L + M
    R = (
      - frac(j+k-i, i) * I2**2/I0**2
      + frac(j+k-1, n) * I4/I0
      - frac(n-1, n) * I22/I0
      )
    H = frac(1, i) * (
        + (j+k-i) * I2**2/I0**2
        + frac(j+k, n) * ((i-1) * I4/I0 - (n-1)*I22/I0)
        )
    Q = L/M*R + H

    G = frac(1, N)
    a = sqrt(frac(n, i) * I2/I0)
    b = sqrt(frac(n, j+k) * (I2/I0 + t * sqrt(k/j*Q)))
    c = sqrt(frac(n, j+k) * (I2/I0 - t * sqrt(j/k*Q)))
    return G, a, b, c
