# -*- coding: utf-8 -*-
#
import numpy
from sympy import Rational as fr, factorial as fact

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
                fr((-1)**i * 2**(-2*s) * (d+n-2*i)**d, fact(i) * fact(d+n-i)),
                numpy.array([
                    [fr(2*p + 1, d+n-2*i) for p in part]
                    for part in partition(s-i, n+1)
                    ])
            )
            for i in range(s+1)
            ]

        self.bary, self.weights = untangle(data)
        self.weights /= sum(self.weights)
        self.points = self.bary[:, 1:]
        return
