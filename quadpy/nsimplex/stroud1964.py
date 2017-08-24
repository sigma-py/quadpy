# -*- coding: utf-8 -*-
#
from __future__ import division

import warnings

import numpy

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
        self.degree = 0

        # TODO find out what's wrong
        warnings.warn('Stroud\'s scheme has degree 0, not 3.')

        # Stroud's book only gives numerical values for certain n, the article
        # explains it in more detail.
        # r is a root of a polynomial of degree 3.
        p = numpy.array([
            n+1, -3, 3/(n+2), -1/(n+2)/(n+3)
            ])
        roots = numpy.sort(numpy.roots(p))
        # all roots are real-valued
        if n > 8:
            assert variant == 'b', 'Choose variant b for n >= 9.'

        r = roots[0] if variant == 'a' else roots[1]

        # s and t are zeros of a polynomial of degree 2
        p = numpy.array([
            1, -(1-(n-1)*r), n/2/(n+2) - (n-1)*r + n*(n-1)/2*r**2
            ])
        s, t = numpy.sort(numpy.roots(p))

        data = [
            (1/n/(n+1), rd(n+1, [(r, 1), (s, 1), (t, n-1)])),
            ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return
