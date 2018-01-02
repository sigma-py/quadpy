# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, fsd, z


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, n, index, symbolic=True):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.dim = n

        if index == '1-n':
            self.degree = 3
            data = [
                (frac(1, 2*n), fsd(n, (sqrt(frac(n, 3)), 1))),
                ]
        else:
            assert index == '2-n'
            self.degree = 5
            r = sqrt(frac(3, 5))
            data = [
                (frac(25*n**2 - 115*n + 162, 162), z(n)),
                (frac(70 - 25*n, 162), fsd(n, (r, 1))),
                (frac(25, 324), fsd(n, (r, 2))),
                ]

        self.points, self.weights = untangle(data)
        reference_volume = 2**n
        self.weights *= reference_volume
        return
