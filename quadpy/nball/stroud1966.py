# -*- coding: utf-8 -*-
#
import numpy

from .helpers import volume_unit_ball

from ..helpers import untangle, fsd, fsd2, pm, combine, z, pm_array


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
        self.name = 'Stroud66'
        self.degree = 5

        if variant == 'a':
            a = numpy.sqrt(2.0 * (n+4))

            r2 = (n + 4 - a)/(n+4)
            s2 = (n*(n+4) + 2*a)/(n**2 + 2*n - 4)/(n+4)

            B1 = 1.0 / (n+2) / (n+4) / r2**2
            B2 = 1.0 / 2**n / (n+2) / (n+4) / s2**2

            r = numpy.sqrt(r2)
            s = numpy.sqrt(s2)

            data = [
                (B1, fsd(n, r, 1)),
                (B2, pm(n, s)),
                ]
        elif variant == 'b':
            alpha = 0.0
            s = numpy.sqrt((n+alpha+2) / (n+2) / (n+alpha+4))
            data = []
            B0 = 1.0
            for k in range(1, n+1):
                B = 2.0**(k-n) * (n+2) * (n+alpha) * (n+alpha+4) \
                    / n / (k+1) / (k+2) / (n+alpha+2)**2
                B0 -= 2.0**(n-k+1) * B
                r = numpy.sqrt((k+2) * (n+alpha+2) / (n+2) / (n+alpha+4))
                v = numpy.concatenate([
                    numpy.zeros((2**(n-k+1), k-1)),
                    pm_array(numpy.array([r] + (n-k) * [s]))
                    ], axis=-1)
                data.append((B, v))
            data.append((B0, z(n)))
        elif variant == 'c':
            a = numpy.sqrt(2.0 * (n+2))
            r = numpy.sqrt((n + 2 + (n-1)*a) / n / (n+4))
            s = numpy.sqrt((n + 2 - a) / n / (n+4))

            B0 = 4.0 / (n+2)**2
            B1 = (n+4) / 2.0**n / (n+2)**2

            data = [
                (B0, z(n)),
                (B1, combine([[+r, -r]] + [[+s, -s]] * (n-1))),
                ]
        elif variant == 'd':
            a = 2*numpy.sqrt(n + 4.0)
            b = numpy.sqrt(2.0 * (n+1) * (n+2) * (n+4))

            r = numpy.sqrt((n*(n+4) + a + (n-1)*b) / n / (n+2) / (n+4))
            s = numpy.sqrt((n*(n+4) + a - b) / n / (n+2) / (n+4))
            t = numpy.sqrt((n + 4 - a) / (n+2) / (n+4))

            B = 1.0 / 2**n / (n+1)

            data = [
                (B, combine([[+r, -r]] + [[+s, -s]] * (n-1))),
                (B, pm(n, t)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= volume_unit_ball(n)
        return
