# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle, fsd, z
from .helpers import volume_unit_ball


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, n, index, alpha):
        self.dim = n

        if index == '11-n':
            self.degree = 3
            r = numpy.sqrt((n + alpha) / (n + alpha + 2))
            data = [
                (1.0/(2*n), fsd(n, r, 1))
                ]
        else:
            assert index == '12-n'
            self.degree = 5
            r = numpy.sqrt(3*(n+alpha+2) / (n+2) / (n+alpha+4))
            B1 = (4.0-n) * (n+2) * (n+alpha) * (n+alpha+4) \
                / (18.0 * n * (n+alpha+2)**2)
            B2 = (n+2.0) * (n+alpha) * (n+alpha+4) \
                / (36.0 * n * (n+alpha+2)**2)
            B0 = 1.0 - 2*n*B1 - 2*n*(n-1)*B2

            data = [
                (B0, z(n)),
                (B1, fsd(n, r, 1)),
                (B2, fsd(n, r, 2)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= volume_unit_ball(n)
        return
