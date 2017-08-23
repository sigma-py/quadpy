# -*- coding: utf-8 -*-
#
from __future__ import division

import math

import numpy

from ..helpers import untangle, rd


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Integration Over Simplexes,
    Mathematical Tables and Other Aids to Computation,
    Vol. 10, No. 55 (Jul., 1956), pp. 137-139,
    <https://doi.org/10.2307/2002484>.
    '''
    def __init__(self, n, index):
        self.dim = n
        if index == '1a':
            self.degree = 2
            r = (n + 2 - math.sqrt(n+2.0)) / (n+1) / (n+2)
            s = (n + 2 + n*math.sqrt(n+2.0)) / (n+1) / (n+2)
            data = [
                (1.0/(n+1), rd(n+1, [(r, n), (s, 1)]))
                ]
        elif index == '1b':
            self.degree = 2
            r = (n + 2 + math.sqrt(n+2.0)) / (n+1) / (n+2)
            s = (n + 2 - n*math.sqrt(n+2.0)) / (n+1) / (n+2)
            data = [
                (1.0/(n+1), rd(n+1, [(r, n), (s, 1)]))
                ]
        else:
            assert index == '2'
            self.degree = 3

            B = -(n+1)**2 / 4 / (n+2)
            C = (n+3)**2 / 4 / (n+1) / (n+2)

            r = 1 / (n+1)
            s = 1 / (n+3)
            t = 3 / (n+3)

            data = [
                (B, numpy.array([numpy.full(n+1, r)])),
                (C, rd(n+1, [(t, 1), (s, n)])),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
