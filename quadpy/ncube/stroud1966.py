# -*- coding: utf-8 -*-
#
import numpy
from sympy import sqrt, Rational as fr, S

from .helpers import _fs11
from ..helpers import untangle, fsd, z, pm


class Stroud1966(object):
    '''
    A.H. Stroud,
    Some Fifth Degree Integration Formulas for Symmetric Regions,
    Mathematics of Computation,
    Vol. 20, No. 93 (Jan., 1966), pp. 90-97,
    Published by: American Mathematical Society,
    <https://doi.org/10.1090/S0025-5718-1966-0191094-8>.
    '''
    def __init__(self, n, variant):
        self.name = 'Stroud66{}'.format(variant)
        self.degree = 5

        if variant == 'a':
            r = sqrt(fr(5*n + 4, 30))
            s = sqrt(fr(5*n + 4, 15*n - 12))
            data = [
                (fr(40, (5*n+4)**2), fsd(n, (r, 1))),
                (fr(5*n - 4, (5*n + 4))**2 / 2**n, pm(n, s)),
                ]
        elif variant == 'b':
            s = 1 / sqrt(3)
            data = [(fr(4, 5*n + 4), z(n))]
            for k in range(1, n+1):
                r = sqrt(fr(5*k + 4, 15))
                arr = numpy.full((2**(n-k+1), n), S.Zero)
                arr[:, k-1:] = pm(n-k+1, 1)
                arr[:, k-1] *= r
                arr[:, k:] *= s
                b = fr(5 * 2**(k-n+1), (5*k-1) * (5*k+4))
                data.append((b, arr))
        elif variant == 'c':
            r = sqrt((5*n + 4 + 2*(n-1)*sqrt(5*n+4)) / (15*n))
            s = sqrt((5*n + 4 - 2*sqrt(5*n+4)) / (15*n))
            data = [
                (fr(4, 5*n+4), z(n)),
                (fr(5, (5*n+4) * 2**n), _fs11(n, r, s)),
                ]
        else:
            assert variant == 'd'
            assert n >= 3
            r = sqrt((5*n - 2*sqrt(5) + 2*(n-1)*sqrt(5*n+5)) / (15*n))
            # This sqrt() is imaginary for negative for n=2.
            s = sqrt((5*n - 2*sqrt(5) - 2*sqrt(5*n+5)) / (15*n))
            t = sqrt((5 + 2*sqrt(5)) / 15)
            w = fr(1, 2**n * (n+1))
            data = [
                (w, _fs11(n, r, s)),
                (w, pm(n, t)),
                ]

        self.points, self.weights = untangle(data)
        reference_volume = 2**n
        self.weights *= reference_volume
        return
