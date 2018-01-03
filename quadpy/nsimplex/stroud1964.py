# -*- coding: utf-8 -*-
#
from __future__ import division

from mpmath import mp
import numpy
import sympy

from ..helpers import untangle, rd


class Stroud1964(object):
    '''
    A.H. Stroud,
    Approximate Integration Formulas of Degree 3 for Simplexes,
    Mathematics of Computation, Vol. 18, No. 88 (Oct., 1964), pp. 590-597,
    <https://doi.org/10.2307/2002945>.
    '''
    def __init__(self, n, variant, symbolic=False):
        assert variant in ['a', 'b']

        frac = sympy.Rational if symbolic else lambda x, y: x/y
        roots = mp.polyroots if symbolic else numpy.roots

        self.dim = n
        self.degree = 3

        if n == 2:
            # The roots sum up to 1
            r, s, t = mp.polyroots([1, -1, frac(1, 4), -frac(1, 60)])
            data = [
                (frac(1, n*(n+1)), rd(n+1, [(r, 1), (s, 1), (t, 1)])),
                ]
        else:
            assert n > 2

            # Stroud's book only gives numerical values for certain n; the
            # article explains it in more detail, namely: r is a root of a
            # polynomial of degree 3.
            rts = numpy.sort(roots([
                n+1, -3, frac(3, n+2), -frac(1, (n+2)*(n+3))
                ]))

            # all roots are real-valued
            if n > 8:
                assert variant == 'b', 'Choose variant b for n >= 9.'

            r = rts[0] if variant == 'a' else rts[1]

            # s and t are zeros of a polynomial of degree 2
            s, t = numpy.sort(roots([
                1,
                -(1-(n-1)*r),
                frac(n, 2*(n+2)) - (n-1)*r + frac(n*(n-1), 2)*r**2
                ]))

            data = [
                (frac(1, n*(n+1)), rd(n+1, [(r, n-1), (s, 1), (t, 1)])),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
