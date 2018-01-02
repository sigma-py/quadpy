# -*- coding: utf-8 -*-
#
from __future__ import division
import numpy
import scipy.special
import sympy

from ..helpers import untangle, fsd, z


class Dobrodeev1970(object):
    '''
    L.N. Dobrodeev,
    Cubature formulas of the seventh order of accuracy for a hypersphere and a
    hypercube,
    USSR Computational Mathematics and Mathematical Physics,
    Volume 10, Issue 1, 1970, Pages 252–253,
    <https://doi.org/10.1016/0041-5553(70)90084-4>.
    '''
    def __init__(self, n, symbolic=True):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        pi = sympy.pi if symbolic else numpy.pi
        sqrt = sympy.sqrt if symbolic else numpy.sqrt
        gamma = sympy.gamma if symbolic else scipy.special.gamma

        self.name = 'Dobrodeev'
        self.degree = 7
        self.dim = n

        A = frac(1, 8)
        B = frac(5-n, 4)
        C = frac((6 - n) * (1 - n**2) + 36, 4*(n + 3))
        D = frac(81, (n + 3) * (n + 6)**2)
        E = frac(45*n**2 + 324*n + 216, n**2 + 12*n + 36) \
            - frac(n * (n**2 - 12*n + 65), 6)

        r = sqrt(frac(3, n+6))
        data = [
            (A, fsd(n, (r, 3))),
            (B, fsd(n, (r, 2))),
            (C, fsd(n, (r, 1))),
            (D, fsd(n, (1, 1))),
            (E, z(n)),
            ]

        self.points, self.weights = untangle(data)

        self.weights /= (
            frac(n, 2) * gamma(frac(n, 2)) / sqrt(pi)**n
            * frac(27 * (n+2) * (n+4), (n+6)**2)
            )
        return
