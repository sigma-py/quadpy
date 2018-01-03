# -*- coding: utf-8 -*-
#
from __future__ import division

from mpmath import mp
import numpy
import sympy


from ..helpers import untangle, rd


class Stroud1966(object):
    '''
    A.H. Stroud,
    Some approximate integration formulas of degree 3 for ann-dimensional
    simplex,
    Numerische Mathematik, November 1966, Volume 9, Issue 1, pp 38–45,
    <https://doi.org/10.1007/BF02165227>.
    '''
    def __init__(self, n, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
        roots = mp.polyroots if symbolic else numpy.roots

        self.dim = n
        self.degree = 3
        if index == 'I':
            sqrt4n13 = sqrt(4*n + 13)

            r = (2*n + 5 - sqrt4n13) / 2 / (n+1) / (n+3)
            s = 1 - n*r

            B = (1 - sqrt4n13) / 2 / (n+1) / (n+2) / (n+3)
            C = (2*n**2 + 10*n + 11 + sqrt4n13) / 2 / (n+1) / (n+2) / (n+3)

            data = [
                (B, rd(n+1, [(1, 1)])),
                (C, rd(n+1, [(r, n), (s, 1)])),
                ]
        elif index == 'II':
            # r is a smallest real-valued root of a polynomial of degree 3
            rts = roots([
                2*(n-2)*(n+1)*(n+3), -(5*n**2+5*n-18), 4*n, - 1
                ])
            r = numpy.min([r.real for r in rts if abs(r.imag) < 1.0e-15])

            s = 1 - n*r
            t = frac(1, 2)

            B = (n-2) / (1 - 2*n*r**2 - 2*(1-n*r)**2) / (n+1) / (n+2)
            C = 2 * (frac(1, n+1) - B) / n

            data = [
                (B, rd(n+1, [(r, n), (s, 1)])),
                (C, rd(n+1, [(t, 2)])),
                ]
        elif index == 'III':
            assert n > 2

            r = frac(1, 2)
            s = frac(1, n)

            prod = (n+1) * (n+2) * (n+3)
            B = frac(6-n, prod)
            C = frac(8*(n-3), (n-2) * prod)
            D = frac(n**3, (n-2) * prod)

            data = [
                (B, rd(n+1, [(1, 1)])),
                (C, rd(n+1, [(r, 2)])),
                (D, rd(n+1, [(s, n)])),
                ]
        elif index == 'IV':
            assert n >= 3

            r = frac(1, n+1)
            s = frac(1, 3)

            A = frac((n+1)**2 * (n-3), (n-2) * (n+2) * (n+3))
            B = frac((9-n), 2 * (n+1) * (n+2) * (n+3))
            C = frac(27, (n-2) * (n+1) * (n+2) * (n+3))

            data = [
                (A, numpy.full((1, n+1), r)),
                (B, rd(n+1, [(1, 1)])),
                (C, rd(n+1, [(s, 3)])),
                ]
        elif index == 'V':
            r = frac(1, n)
            s = frac(1, 3)

            prod = (n+1) * (n+2) * (n+3)
            A = frac(-n**2 + 11*n - 12, 2 * (n-1) * prod)
            B = frac(n**3, (n-1) * prod)
            C = frac(27, (n-1) * prod)

            data = [
                (A, rd(n+1, [(1, 1)])),
                (B, rd(n+1, [(r, n)])),
                (C, rd(n+1, [(s, 3)])),
                ]
        elif index == 'VI':
            assert n >= 3
            assert n != 5

            r = frac(1, n+1)
            s = frac(1, 3)
            t = frac(1, n-2)

            prod = (n+1) * (n+2) * (n+3)
            A = frac((3-n) * (n-12) * (n+1)**2, 3 * (n-2) * (n+2) * (n+3))
            B = frac(54 * (3*n-11), (n-5) * (n-2) * (n-1) * prod)
            C = frac(2 * (n-2)**2 * (n-9), (n-5) * (n-1) * prod)

            data = [
                (A, numpy.full((1, n+1), r)),
                (B, rd(n+1, [(s, 3)])),
                (C, rd(n+1, [(t, n-2)])),
                ]
        else:
            assert index == 'VII'
            assert n >= 3
            assert n != 5

            s = frac(1, 3)
            t = frac(1, n-2)

            prod = (n+1) * (n+2) * (n+3)
            A = frac((12-n), 2 * prod)
            B = frac(27 * (n-7), (n-5) * (n-1) * prod)
            C = frac(6 * (n-2)**2, (n-5) * (n-1) * prod)

            data = [
                (A, rd(n+1, [(1, 1)])),
                (B, rd(n+1, [(s, 3)])),
                (C, rd(n+1, [(t, n-2)])),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
