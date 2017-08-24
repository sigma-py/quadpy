# -*- coding: utf-8 -*-
#
from __future__ import division

import orthopy
import numpy

from ..helpers import untangle
from ..line_segment.gauss_legendre import GaussLegendre
from ..nsimplex.stroud import Stroud as nsimplex_Stroud


class Stroud(object):
    '''
    A. H. Stroud,
    Approximate calculation of multiple integrals,
    Englewood Cliffs, NJ : Prentice-Hall, c 1971. - XIII,
    ISBN 0-13-043893-6.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        self.name = 'Stroud({})'.format(index)
        if index == 'T3 5-1':
            self.set_data(nsimplex_Stroud(3, 'Tn 5-1'))
        else:
            assert index == 'T3 7-1'
            self.degree = 7

            gl4 = GaussLegendre(4)
            r = (gl4.points + 1) / 2
            A = gl4.weights / 2

            # Generate Gauss formula for int_0^1 (1-s) * f(s) ds.
            k = numpy.arange(8)
            moments = 1 / (k**2 + 3*k + 2)
            alpha, beta = orthopy.chebyshev(moments)
            s, B = orthopy.gauss_from_coefficients(alpha, beta)

            # Generate Gauss formula for int_0^1 (1-t)^2 * f(t) ds.
            k = numpy.arange(8)
            moments = 2 / (k**3 + 6*k**2 + 11*k + 6)
            alpha, beta = orthopy.chebyshev(moments)
            t, C = orthopy.gauss_from_coefficients(alpha, beta)

            data = [
                (6*A[i]*B[j]*C[k], numpy.array([[
                    t[k],
                    s[j]*(1-t[k]),
                    r[i]*(1-s[j])*(1-t[k]),
                    (1-r[i])*(1-s[j])*(1-t[k])
                    ]]))
                for i in range(4)
                for j in range(4)
                for k in range(4)
                ]

            self.bary, self.weights = untangle(data)
            self.points = self.bary[:, 1:]
        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        self.bary = scheme.bary
        return
