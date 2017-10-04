# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr, binomial
import warnings

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
    def __init__(self, n):
        warnings.warn('The Phillips schemes are only single-precision.')
        self.name = 'Phillips'
        self.degree = 7

        if n == 2:
            lambda1 = 1.0
            lambda2 = 0.462910050
            mu = 0.774596669
            a = -0.158024691
            b1 = 0.036363636
            b2 = 0.175982043
            c = 0.077160494
        elif n == 3:
            lambda1 = 1.0
            lambda2 = 0.597614305
            mu = 0.632455532
            nu = 1.0
            a = 0.542962963
            b1 = 0.032098765
            b2 = -0.193580247
            c = 0.115740741
            d = 0.004629630
        elif n == 4:
            lambda1 = 1.0
            lambda2 = 0.313391585
            mu = 0.447213596
            nu = 0.707106781
            a = -10.888215488
            b1 = 0.025082508
            b2 = 2.007240724
            c = -0.231481481
            d = 0.037037037
        else:
            assert n >= 5
            p1 = 1
            En = fr(25*n**2 - 165*n + 302, 972)
            p2 = 1 / (fr(3, 5) - fr(1, 35*En))

            a1 = fr(3, 5) * En
            a2 = fr(9, 25) * En + fr(2, 175)
            beta1 = (a1 - a2 * p2) / (p1 - p2)
            beta2 = (a1 - a2 * p1) / (p2 - p1)

            lambda1 = 1 / sqrt(p1)
            lambda2 = 1 / sqrt(p2)
            mu = sqrt(fr(3, 5))
            nu = sqrt(fr(3, 5))

            b1 = beta1 / lambda1**6
            b2 = beta2 / lambda2**6

            gamma = fr((n-1) * (19 - 5*n), 270)
            delta = fr((n-1) * (n - 2), 108)

            c = gamma / (2 * (n-1) * mu**6)
            d = delta / (4 * binomial(n-1, 2) * nu**6)

            a = 1 - 2*n*(b1+b2) - 4*binomial(n, 2) * c - 8*binomial(n, 3) * d

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
