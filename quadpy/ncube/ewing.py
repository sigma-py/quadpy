# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from ..helpers import untangle, z, pm


class Ewing(object):
    '''
    G.M. Ewing,
    On Approximate Cubature,
    The American Mathematical Monthly,
    Vol. 48, No. 2 (Feb., 1941), pp. 134-136,
    <https://doi.org/10.2307/2303604>.
    '''
    def __init__(self, n, symbolic=True):
        frac = sympy.Rational if symbolic else lambda x, y: x/y

        self.degree = 3
        data = [
            (frac(2, 3), z(n)),
            (frac(1, 3 * 2**n), pm(n, 1)),
            ]

        self.points, self.weights = untangle(data)
        reference_volume = 2**n
        self.weights *= reference_volume
        return
