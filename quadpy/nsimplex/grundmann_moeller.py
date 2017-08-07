# -*- coding: utf-8 -*-
#
from math import factorial
import numpy

from ..helpers import partition, untangle


class GrundmannMoeller(object):
    '''
    A. Grundmann and H.M. Moeller,
    Invariant integration formulas for the n-simplex by combinatorial methods,
    SIAM J. Numer. Anal. 15 (1978), 282-290,
    <http://dx.doi.org/10.1137/0715019>.

    Abstract:
    For the n-simplex T_n, integration formulas of arbitrary odd degree are
    derived and the monomial representations of the orthogonal polynomials
    corresponding to T_n are given.
    '''
    def __init__(self, n, s):
        self.name = 'GrundmannMöller(dim={}, {})'.format(n, s)
        d = 2*s + 1
        self.degree = d
        self.dim = n

        data = [
            (
                (-1)**i * 2.0**(-2.0*s) * (d + n - 2*i)**d
                / factorial(i) / factorial(d + n - i),
                (2.0*numpy.array(partition(s-i, n+1)) + 1.0) / (d+n-2*i)
            )
            for i in range(s+1)
            ]

        self.bary, self.weights = untangle(data)
        self.weights /= numpy.sum(self.weights)
        self.points = self.bary[:, 1:]
        return
