# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from ..helpers import untangle, rd


class Stroud1961(object):
    '''
    A.H. Stroud,
    Numerical Integration Formulas of Degree 3 for Product Regions and Cones
    Mathematics of Computation, Vol. 15, No. 74 (Apr., 1961), pp. 143-150,
    <https://doi.org/10.2307/2004220>.
    '''
    def __init__(self, n):
        self.dim = n
        self.degree = 3

        r = 1 / (n+1)
        s = 1 / n

        A = (3-n) * (n+1)**2 / (n+2) / (n+3)
        B = 3 / (n+1) / (n+2) / (n+3)
        C = n**3 / (n+1) / (n+2) / (n+3)

        data = [
            (A, numpy.array([numpy.full(n+1, r)])),
            (B, rd(n+1, [(1.0, 1)])),
            (C, rd(n+1, [(s, n)])),
            ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
