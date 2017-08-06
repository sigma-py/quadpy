# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _fsd, _z
from ..helpers import untangle


class Dobrodeev(object):
    '''
    L.N. Dobrodeev,
    Cubature formulas of the seventh order of accuracy for a hypersphere and a
    hypercube,
    USSR Computational Mathematics and Mathematical Physics,
    Volume 10, Issue 1, 1970, Pages 252–253,
    <https://doi.org/10.1016/0041-5553(70)90084-4>.
    '''
    def __init__(self, n):
        self.name = 'Dobrodeev'
        self.degree = 7
        self.dim = n

        A = 0.125
        B = (5.0 - n) / 4.0
        C = ((6.0 - n) * (1.0 - n**2) + 36.0) / 4.0 / (n + 3.0)
        D = 81.0 / (n + 3.0) / (n + 6.0)**2
        E = (45.0*n**2 + 324.0*n + 216.0) / (n**2 + 12.0*n + 36.0) \
            - n * (n**2 - 12.0*n + 65.0) / 6.0

        r = numpy.sqrt(3.0 / (n + 6.0))
        data = [
            (A, _fsd(n, r, 3)),
            (B, _fsd(n, r, 2)),
            (C, _fsd(n, r, 1)),
            (D, _fsd(n, 1.0, 1)),
            (E, _z(n)),
            ]

        self.points, self.weights = untangle(data)

        self.weights /= (
            (0.5*n) * numpy.math.gamma(0.5*n) / numpy.pi**(0.5*n)
            * 27.0 * (n+2.0) * (n+4.0) / (n+6.0)**2
            )
        return
