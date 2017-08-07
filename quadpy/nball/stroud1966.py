# -*- coding: utf-8 -*-
#
import numpy

from .helpers import volume_unit_ball

from ..helpers import untangle, fsd, pm


class Stroud1966(object):
    '''
    A.H. Stroud,
    Some Fifth Degree Integration Formulas for Symmetric Regions,
    Mathematics of Computation,
    Vol. 20, No. 93 (Jan., 1966), pp. 90-97,
    Published by: American Mathematical Society,
    <https://doi.org/10.1090/S0025-5718-1966-0191094-8>.
    '''
    def __init__(self, n):
        self.name = 'Stroud66'
        self.degree = 5

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

        self.points, self.weights = untangle(data)
        print(sum(self.weights))
        self.weights *= volume_unit_ball(n)
        return
