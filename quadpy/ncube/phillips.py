# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import scipy.special
import sympy

from ..helpers import untangle, fsd, z


class Phillips(object):
    '''
    G.M. Phillips,
    Numerical integration over an N-dimensional rectangular region,
    Comput J (1967) 10 (3): 297-299,
    <https://doi.org/10.1093/comjnl/10.3.297>.

    Abstract:
    Gaussian-type formulae are derived for all values of N >= 2.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = 'Phillips'
        self.degree = 7

        if n == 2:
            p1 = 1
            p2 = frac(14, 3)
            q = frac(5, 3)
        elif n == 3:
            p1 = 1
            p2 = frac(14, 5)
            q = frac(5, 2)
            r = 1
        elif n == 4:
            p1 = 1
            p2 = frac(112, 11)
            q = 5
            r = 2
        else:
            assert n >= 5
            p1 = 1
            En = frac(25*n**2 - 165*n + 302, 972)
            p2 = 1 / (frac(3, 5) - frac(1, 35*En))
            q = frac(5, 3)
            r = frac(5, 3)

        gamma = frac((n-1) * (19 - 5*n), 270)
        delta = frac((n-1) * (n - 2), 108)

        a1 = frac(23 - 5*n, 180) - gamma*q / 2
        a2 = frac(35*n**2 - 231*n + 466, 3780)
        beta1 = (a1 - a2 * p2) / (p1 - p2)
        beta2 = (a1 - a2 * p1) / (p2 - p1)

        lambda1 = 1 / sqrt(p1)
        lambda2 = 1 / sqrt(p2)
        mu = 1 / sqrt(q)

        b1 = beta1 / lambda1**6
        b2 = beta2 / lambda2**6

        c = gamma / (2 * (n-1) * mu**6)

        a = 1 - 2*n*(b1+b2) - 4*int(scipy.special.binom(n, 2)) * c

        if n > 2:
            nu = 1 / sqrt(r)
            d = delta / (4 * int(scipy.special.binom(n-1, 2)) * nu**6)
            a -= 8*int(scipy.special.binom(n, 3)) * d

        data = [
            (a, z(n)),
            (b1, fsd(n, (lambda1, 1))),
            (b2, fsd(n, (lambda2, 1))),
            (c, fsd(n, (mu, 2))),
            ]

        if n > 2:
            data.append((d, fsd(n, (nu, 3))))

        self.points, self.weights = untangle(data)
        reference_volume = 2**n
        self.weights *= reference_volume
        return
