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
        self.degree = 3
        if index == 'I':
            sqrt4n13 = sqrt(4*n+13)

            r = (2*n + 5 - sqrt4n13) / 2 / (n+1) / (n+3)
            s = 1 - n*r

            B = (1 - sqrt4n13) / 2 / (n+1) / (n+2) / (n+3)
            C = (2*n**2 + 10*n + 11 + sqrt4n13) / 2 / (n+1) / (n+2) / (n+3)

            data = [
                (B, rd(n+1, [(1.0, 1)])),
                (C, rd(n+1, [(r, n), (s, 1)])),
                ]
        elif index == 'II':
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
        elif index == 'III':
            assert n > 2

            r = 0.5
            s = 1 / n

            B = (6-n) / (n+1) / (n+2) / (n+3)
            C = 8*(n-3) / (n-2) / (n+1) / (n+2) / (n+3)
            D = n**3 / (n-2) / (n+1) / (n+2) / (n+3)

            data = [
                (B, rd(n+1, [(1.0, 1)])),
                (C, rd(n+1, [(r, 2)])),
                (D, rd(n+1, [(s, n)])),
                ]
        elif index == 'IV':
            assert n >= 3

            r = 1/(n+1)
            s = 1/3

            A = (n+1)**2 * (n-3) / (n-2) / (n+2) / (n+3)
            B = (9-n) / 2 / (n+1) / (n+2) / (n+3)
            C = 27 / (n-2) / (n+1) / (n+2) / (n+3)

            data = [
                (A, numpy.array([numpy.full(n+1, r)])),
                (B, rd(n+1, [(1.0, 1)])),
                (C, rd(n+1, [(s, 3)])),
                ]
        elif index == 'V':
            r = 1/n
            s = 1/3

            A = (-n**2 + 11*n - 12) / 2 / (n-1) / (n+1) / (n+2) / (n+3)
            B = n**3 / (n-1) / (n+1) / (n+2) / (n+3)
            C = 27 / (n-1) / (n+1) / (n+2) / (n+3)

            data = [
                (A, rd(n+1, [(1.0, 1)])),
                (B, rd(n+1, [(r, n)])),
                (C, rd(n+1, [(s, 3)])),
                ]
        elif index == 'VI':
            assert n >= 3
            assert n != 5

            r = 1/(n+1)
            s = 1/3
            t = 1/(n-2)

            A = (3-n) * (n-12) * (n+1)**2 / 3 / (n-2) / (n+2) / (n+3)
            B = 54 * (3*n-11) / (n-5) / (n-2) / (n-1) / (n+1) / (n+2) / (n+3)
            C = 2 * (n-2)**2 * (n-9) / (n-5) / (n-1) / (n+1) / (n+2) / (n+3)

            data = [
                (A, numpy.array([numpy.full(n+1, r)])),
                (B, rd(n+1, [(s, 3)])),
                (C, rd(n+1, [(t, n-2)])),
                ]
        else:
            assert index == 'VII'
            assert n >= 3
            assert n != 5

            s = 1/3
            t = 1/(n-2)

            A = (12-n) / 2 / (n+1) / (n+2) / (n+3)
            B = 27 * (n-7) / (n-5) / (n-1) / (n+1) / (n+2) / (n+3)
            C = 6 * (n-2) * (n-2) / (n-5) / (n-1) / (n+1) / (n+2) / (n+3)

            data = [
                (A, rd(n+1, [(1.0, 1)])),
                (B, rd(n+1, [(s, 3)])),
                (C, rd(n+1, [(t, n-2)])),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
