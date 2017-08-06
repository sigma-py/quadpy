# -*- coding: utf-8 -*-
#
from __future__ import division

import math
from math import factorial as fact, sqrt
import numpy

from ..helpers import untangle, fsd, fsd2


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
    def __init__(self, n):
        self.name = 'Dobrodeev1978'
        self.degree = 5
        self.dim = n

        dim_config = {
            2: ('II', 1, 1, 1),
            3: ('II', 1, 1, 1),
            4: ('II', 4, 1, 1),
            5: ('I', 3, 1, 1),
            6: ('I', 3, 1, 1),
            7: ('I', 7, 1, 1),
            8: ('I', 8, 1, 1),
            9: ('I', 9, 1, 1),
            10: ('I', 9, 1, 2),
            11: ('I', 5, 1, 2),
            12: ('I', 12, 1, 2),
            13: ('I', 13, 1, 2),
            14: ('I', 14, 1, 2),
            15: ('I', 15, 1, 2),
            16: ('II', 15, 1, 4),
            17: ('I', 17, 1, 3),
            18: ('I', 18, 1, 3),
            19: ('I', 19, 1, 3),
            20: ('I', 20, 1, 3),
            }

        pm_type, i, j, k = dim_config[n]

        t = 1 if pm_type == 'I' else -1

        I0 = _integrate_monomial_over_unit_nball(n * [0])
        I2 = _integrate_monomial_over_unit_nball([2] + (n-1) * [0])
        I22 = _integrate_monomial_over_unit_nball([2, 2] + (n-2) * [0])
        I4 = _integrate_monomial_over_unit_nball([4] + (n-1) * [0])

        L = fact(n) // (fact(i) * fact(n-i)) * 2**i
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

        data = [
            (G, fsd(n, a, i)),
            (G, fsd2(n, b, c, j, k)),
            ]

        self.points, self.weights = untangle(data)

        self.weights *= I0
        return


def _integrate_monomial_over_unit_nball(exp):
    '''
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448,
    <https://dx.doi.org/10.2307/2695802>.
    '''
    radius = 1.0
    n = len(exp)
    alpha = n + sum(exp)
    return radius**alpha / alpha * _integrate_monomial_over_unit_sphere(exp)


def _integrate_monomial_over_unit_sphere(alpha):
    '''
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448,
    <https://dx.doi.org/10.2307/2695802>.
    '''
    alpha = numpy.array(alpha)
    if any(alpha % 2 == 1):
        return 0.0
    # Use lgamma since other with ordinary gamma, numerator and denominator
    # might overflow.
    return 2.0 * math.exp(
        math.fsum([math.lgamma(0.5*(a+1)) for a in alpha])
        - math.lgamma(math.fsum([0.5*(a+1) for a in alpha]))
        )
