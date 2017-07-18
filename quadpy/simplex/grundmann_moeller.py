# -*- coding: utf-8 -*-
#
from math import factorial
import numpy


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

        bary = []
        weights = []
        for i in range(s+1):
            beta = numpy.array(list(_weak_compositions(n+1, s-i)))
            bary.append(
                (2.0*beta + 1.0) / (d+n-2*i)
                )
            weights.append(
                numpy.full(
                    len(beta),
                    (-1)**i * 2.0**(-2.0*s) * (d + n - 2*i)**d
                    / factorial(i) / factorial(d + n - i)
                    ))

        self.weights = numpy.concatenate(weights)
        self.weights /= numpy.sum(self.weights)

        self.bary = numpy.concatenate(bary)
        self.points = self.bary[:, 1:]
        return


def _weak_compositions(boxes, balls, parent=tuple()):
    # https://stackoverflow.com/a/36748940/353337
    if boxes > 1:
        for i in range(balls + 1):
            for x in _weak_compositions(boxes - 1, i, parent + (balls - i,)):
                yield x
    else:
        yield parent + (balls,)
