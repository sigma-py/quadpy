# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
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

    def __init__(self, n, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = 'Dobrodeev'
        self.degree = 7
        self.dim = n

        A = frac(1, 8)
        B = frac(19-5*n, 20)
        alpha = 35*n * (5*n - 33)
        C = frac((alpha + 2114)**3, 700 * (alpha+1790.0) * (alpha+2600.0))
        D = frac(729, 1750) * frac(alpha + 2114, alpha + 2600)
        E = frac(n * (n-1) * (n - 4.7), 3) - 2*n * (C + D) + frac(729, 125)

        a = sqrt(frac(3, 5))
        b = a
        c = sqrt(frac(3, 5) * frac(alpha+1790, alpha+2114))
        data = [
            (A, fsd(n, (a, 3))),
            (B, fsd(n, (b, 2))),
            (C, fsd(n, (c, 1))),
            (D, fsd(n, (1.0, 1))),
            (E, z(n)),
            ]

        self.points, self.weights = untangle(data)
        self.weights /= frac(729, 125 * 2**n)
        return
