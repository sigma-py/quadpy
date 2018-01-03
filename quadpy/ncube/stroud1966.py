# -*- coding: utf-8 -*-
#
import numpy
import sympy

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
    # pylint: disable=too-many-locals
    def __init__(self, n, variant, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt
        zero = sympy.S.Zero if symbolic else 0.0

        self.name = 'Stroud66{}'.format(variant)
        self.degree = 5

        if variant == 'a':
            r = sqrt(frac(5*n + 4, 30))
            s = sqrt(frac(5*n + 4, 15*n - 12))
            data = [
                (frac(40, (5*n+4)**2), fsd(n, (r, 1))),
                (frac(5*n - 4, (5*n + 4))**2 / 2**n, pm(n, s)),
                ]
        elif variant == 'b':
            s = 1 / sqrt(3)
            data = [(frac(4, 5*n + 4), z(n))]
            for k in range(1, n+1):
                r = sqrt(frac(5*k + 4, 15))
                arr = numpy.full((2**(n-k+1), n), zero)
                arr[:, k-1:] = pm(n-k+1, 1)
                arr[:, k-1] *= r
                arr[:, k:] *= s
                b = frac(5 * 2**(k-n+1), (5*k-1) * (5*k+4))
                data.append((b, arr))
        elif variant == 'c':
            r = sqrt((5*n + 4 + 2*(n-1)*sqrt(5*n+4)) / (15*n))
            s = sqrt((5*n + 4 - 2*sqrt(5*n+4)) / (15*n))
            data = [
                (frac(4, 5*n+4), z(n)),
                (frac(5, (5*n+4) * 2**n), _fs11(n, r, s)),
                ]
        else:
            assert variant == 'd'
            assert n >= 3
            r = sqrt((5*n - 2*sqrt(5) + 2*(n-1)*sqrt(5*n+5)) / (15*n))
            # This sqrt() is imaginary for negative for n=2.
            s = sqrt((5*n - 2*sqrt(5) - 2*sqrt(5*n+5)) / (15*n))
            t = sqrt((5 + 2*sqrt(5)) / 15)
            w = frac(1, 2**n * (n+1))
            data = [
                (w, _fs11(n, r, s)),
                (w, pm(n, t)),
                ]

        self.points, self.weights = untangle(data)
        reference_volume = 2**n
        self.weights *= reference_volume
        return
