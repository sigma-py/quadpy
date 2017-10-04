# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
from mpmath import mp
from sympy import Rational as fr

from ..helpers import untangle, rd


class Stroud1964(object):
    '''
    A.H. Stroud,
    Approximate Integration Formulas of Degree 3 for Simplexes,
    Mathematics of Computation, Vol. 18, No. 88 (Oct., 1964), pp. 590-597,
    <https://doi.org/10.2307/2002945>.
    '''
    def __init__(self, n, variant):
        assert variant in ['a', 'b']

        self.dim = n
        self.degree = 3

        if n == 2:
            # The roots sum up to 1
            r, s, t = mp.polyroots([1, -1, fr(1, 4), -fr(1, 60)])
            data = [
                (fr(1, n*(n+1)), rd(n+1, [(r, 1), (s, 1), (t, 1)])),
                ]
        else:
            assert n > 2

            # Stroud's book only gives numerical values for certain n, the
            # article explains it in more detail.
            # r is a root of a polynomial of degree 3.
            p = [n+1, -3, fr(3, n+2), -fr(1, (n+2)*(n+3))]

            # mpmath polyroots takes really long, unfortunately.
            # from mpmath import mp
            # roots = mp.polyroots(p)

            roots = numpy.sort(numpy.roots(p))
            # all roots are real-valued
            if n > 8:
                assert variant == 'b', 'Choose variant b for n >= 9.'

            r = roots[0] if variant == 'a' else roots[1]

            # s and t are zeros of a polynomial of degree 2
            # b = 1 - (n-1) * r
            p = [
                1, -(1-(n-1)*r), fr(n, 2*(n+2)) - (n-1)*r + fr(n*(n-1), 2)*r**2
                ]
            s, t = numpy.sort(numpy.roots(p))

            # TODO check what's given wrong in Stroud's book
            data = [
                (fr(1, n*(n+1)), rd(n+1, [(r, n-1), (s, 1), (t, 1)])),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
