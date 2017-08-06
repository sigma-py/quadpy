# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle, fsd, z


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
        B = (19.0 - 5.0*n) / 20.0
        alpha = 35.0*n * (5.0*n - 33.0)
        C = (alpha + 2114.0)**3 / (700.0 * (alpha+1790.0) * (alpha+2600.0))
        D = 729.0/1750.0 * (alpha + 2114.0)/(alpha + 2600.0)
        E = n * (n-1) * (n - 4.7) / 3.0 - 2*n * (C + D) + 729.0/125.0

        a = numpy.sqrt(3.0 / 5.0)
        b = a
        c = numpy.sqrt(3.0 / 5.0 * (alpha+1790.0) / (alpha+2114.0))
        data = [
            (A, fsd(n, a, 3)),
            (B, fsd(n, b, 2)),
            (C, fsd(n, c, 1)),
            (D, fsd(n, 1.0, 1)),
            (E, z(n)),
            ]

        self.points, self.weights = untangle(data)

        self.weights /= 729.0/125.0 / 2**n
        return
