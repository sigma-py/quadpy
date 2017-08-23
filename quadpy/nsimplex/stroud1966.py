# -*- coding: utf-8 -*-
#
from __future__ import division

from math import sqrt

import numpy

from ..helpers import untangle, rd


class Stroud1966(object):
    '''
    A.H. Stroud,
    Some approximate integration formulas of degree 3 for ann-dimensional
    simplex,
    Numerische Mathematik, November 1966, Volume 9, Issue 1, pp 38–45,
    <https://doi.org/10.1007/BF02165227>.
    '''
    def __init__(self, n, index):
        self.dim = n
        if index == 1:
            self.degree = 3

            sqrt4n13 = sqrt(4*n+13)

            r = (2*n + 5 - sqrt4n13) / 2 / (n+1) / (n+3)
            s = 1 - n*r

            B = (1 - sqrt4n13) / 2 / (n+1) / (n+2) / (n+3)
            C = (2*n**2 + 10*n + 11 + sqrt4n13) / 2 / (n+1) / (n+2) / (n+3)

            data = [
                (B, rd(n+1, [(1.0, 1)])),
                (C, rd(n+1, [(r, n), (s, 1)])),
                ]
        elif index == 2:
            self.degree = 3

            # r is a smallest real-valued root of a polynomial of degree 3
            p = numpy.array([
                2*(n-2)*(n+1)*(n+3), -(5*n**2+5*n-18), 4*n, - 1
                ])
            roots = numpy.roots(p)
            i = numpy.where(abs(roots.imag) < 1.0e-12)[0]
            r = numpy.min(numpy.real(roots[i]))

            s = 1 - n*r
            t = 0.5

            B = (n-2) / (n+1) / (n+2) / (1 - 2*n*r**2 - 2*(1-n*r)**2)
            C = 2/n * (1/(n+1) - B)

            data = [
                (B, rd(n+1, [(r, n), (s, 1)])),
                (C, rd(n+1, [(t, 2)])),
                ]
        else:
            assert False

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
