# -*- coding: utf-8 -*-
#
import warnings

import numpy
from sympy import Rational as fr, sqrt

from ..helpers import untangle


class Hillion(object):
    '''
    P. Hillion,
    Numerical Integration on a Triangle,
    International Journal for Numerical Methods in Engineering,
    Vol. 11, 797-815 (1977).
    DOI:10.1002/nme.1620110504,
    <https://dx.doi.org/10.1002/nme.1620110504>.

    Note that the schemes here are not fully symmetric.
    '''
    def __init__(self, index):
        self.name = 'Hillion(%d)' % index

        warnings.warn('Some Hillion schemes are single precision only.')
        # TODO generate more digits

        if index == 1:
            self.degree = 1
            data = [
                (fr(1, 2), [[fr(1, 3), fr(1, 3)]])
                ]
        elif index == 2:
            self.degree = 2
            data = [
                (fr(1, 6), _symm(0, fr(1, 2))),
                (fr(1, 6), numpy.full((1, 2), fr(1, 2))),
                ]
        elif index == 3:
            self.degree = 2
            data = [
                (fr(1, 6), fr(2, 3) - _symm(0, fr(1, 2))),
                (fr(1, 6), fr(2, 3) - numpy.full((1, 2), fr(1, 2))),
                ]
        elif index == 4:
            self.degree = 2
            a0, a1 = [(3 + i*sqrt(3)) / 8 for i in [+1, -1]]
            data = [
                (fr(1, 18), numpy.array([[0, 0]])),
                (fr(2, 9), _symm(a0, a1)),
                ]
        elif index == 5:
            self.degree = 2
            a0, a1 = [(3 + i*sqrt(3)) / 8 for i in [+1, -1]]
            data = [
                (fr(1, 18), fr(2, 3) - numpy.array([[0, 0]])),
                (fr(2, 9), fr(2, 3) - _symm(a0, a1))
                ]
        elif index == 6:
            self.degree = 2
            lmbda = 0.655308609
            mu = 0.247060398
            data = [
                (fr(1, 8), _symm(lmbda, mu)),
                (fr(1, 8), fr(2, 3) - _symm(lmbda, mu)),
                ]
        # Not working
        # elif index == 7:
        #     self.degree = 3
        #     data = [
        #         (0.159020691, _symm(0.666390246, 0.280019915)),
        #         (0.090979309, _symm(0.178558728, 0.075031109)),
        #         ]
        elif index == 8:
            self.degree = 3
            lambda2, lambda3 = [(32 + i*2*sqrt(46))/105 for i in [+1, -1]]
            w1, w2 = [(3266 + i * 19*sqrt(46)) / 17664 for i in [+1, -1]]
            data = [
                (fr(25, 384), _symm(0, fr(4, 5))),
                (w1, numpy.array([[lambda2, lambda2]])),
                (w2, numpy.array([[lambda3, lambda3]])),
                ]
        # Not working. The weights don't even add up.
        # elif index == 9:
        #     self.degree = 3
        #     data = [
        #         (9.0/32.0, numpy.array([[1.0/3.0, 1.0/3.0]])),
        #         (25.0/96.0, _symm(0.2, 0.6)),
        #         (25.0/96.0, numpy.array([[0.2, 0.2]])),
        #         ]
        else:
            assert index == 10
            self.degree = 3
            lambda1, lambda2 = [(16 + i*2*sqrt(14))/25 for i in [+1, -1]]
            w1, w2 = [(161 + i * 17*sqrt(14))/2688 for i in [+1, -1]]
            data = [
                (w2, _symm(lambda1, 0)),
                (w1, _symm(0, lambda2)),
                (fr(25, 96), numpy.array([[fr(2, 5), fr(2, 5)]])),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 2
        return


def _symm(a, b):
    return numpy.array([
        [a, b],
        [b, a],
        ])
