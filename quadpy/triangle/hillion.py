# -*- coding: utf-8 -*-
#
from sympy import Rational as fr, sqrt
import warnings

import numpy

from ..helpers import untangle


class Hillion(object):
    '''
    P. Hillion,
    Numerical Integration on a Triangle,
    International Journal for Numerical Methods in Engineering,
    Vol. 11, 797-815 (1977).
    DOI:10.1002/nme.1620110504,
    <https://dx.doi.org/10.1002/nme.1620110504>.

    Note that the schemes here are not fully symmetric. Also note that in the
    article, the quadrature constants are specified with low precision such
    that the tests are failing. What is needed here is a reimplementation of
    Hillion's method to retrieve more digits.
    '''
    def __init__(self, index):
        self.name = 'Hillion(%d)' % index

        warnings.warn('Some Hillion schemes are single precision only.')

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
            a0 = (3 + sqrt(3)) / 8
            a1 = (3 - sqrt(3)) / 8
            data = [
                (fr(1, 18), numpy.array([[0, 0]])),
                (fr(2, 9), _symm(a0, a1)),
                ]
        elif index == 5:
            self.degree = 2
            a0 = (3 + sqrt(3)) / 8
            a1 = (3 - sqrt(3)) / 8
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
            lambda2 = 0.433949142
            lambda3 = 0.175574667
            data = [
                (0.065104166, _symm(0, fr(4, 5))),
                (0.192191138, numpy.array([[lambda2, lambda2]])),
                (0.177600528, numpy.array([[lambda3, lambda3]])),
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
            data = [
                (0.036232077, _symm(0.939332590, 0)),
                (0.083559589, _symm(0, 0.340667409)),
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
