# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _s21, _s3
from ..helpers import untangle


# pylint: disable=too-many-locals
class Hillion(object):
    '''
    P. Hillion,
    Numerical Integration on a Triangle,
    International Journal for Numerical Methods in Engineering,
    Vol. 11, 797-815 (1977).
    DOI:10.1002/nme.1620110504,
    <https://doi.org/10.1002/nme.1620110504>.

    Note that the schemes here are not fully symmetric.
    '''
    def __init__(self, index, symbolic=False):
        # ENH in the article, most schemes are given only in single precision.
        # quadpy adds symbolic expressions
        frac = sympy.frac if symbolic else lambda x, y: x/y
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        self.name = 'Hillion(%d)' % index

        if index == 1:
            self.degree = 1
            data = [
                (frac(1, 2), _s3(symbolic)),
                ]
        elif index == 2:
            self.degree = 2
            data = [
                (frac(1, 6), _s21(frac(1, 2))),
                ]
        elif index == 3:
            self.degree = 2
            data = [
                (frac(1, 6), _s21(frac(1, 6))),
                ]
        elif index == 4:
            self.degree = 2
            a0, a1 = [(3 + i*sqrt(3)) / 8 for i in [+1, -1]]
            data = [
                (frac(1, 18), numpy.array([[0, 0, 1]])),
                (frac(2, 9), _symm(a0, a1)),
                ]
        elif index == 5:
            self.degree = 2
            a0, a1 = [(3 + i*sqrt(3)) / 8 for i in [+1, -1]]
            data = [
                (
                    frac(1, 18),
                    numpy.array([[frac(2, 3), frac(2, 3), -frac(1, 3)]])
                ),
                (frac(2, 9), _symm(frac(2, 3) - a0, frac(2, 3) - a1))
                ]
        elif index == 6:
            self.degree = 2
            lm, mu = [(2 + i*sqrt(2 + i*sqrt(3))) / 6 for i in [+1, -1]]
            data = [
                (frac(1, 8), _symm(lm, mu)),
                (frac(1, 8), _symm(frac(2, 3) - lm, frac(2, 3) - mu)),
                ]
        elif index == 7:
            self.degree = 3
            a = (6 + sqrt(2) + sqrt(6 * (3 + 2*sqrt(2)))) / 20
            b = (6 + sqrt(2) - sqrt(6 * (3 + 2*sqrt(2)))) / 20
            c = (6 - sqrt(2) + sqrt(6 * (3 - 2*sqrt(2)))) / 20
            d = (6 - sqrt(2) - sqrt(6 * (3 - 2*sqrt(2)))) / 20
            w1 = (2 - 3 * (b+c)) / 12 / (a+d-b-c)
            w2 = (2 - 3 * (a+d)) / 12 / (b+c-a-d)
            data = [
                (w1, _symm(a, d)),
                (w2, _symm(c, b)),
                ]
        elif index == 8:
            self.degree = 3
            lambda2, lambda3 = [(32 + i*2*sqrt(46))/105 for i in [+1, -1]]
            w1, w2 = [(3266 + i * 19*sqrt(46)) / 17664 for i in [+1, -1]]
            data = [
                (frac(25, 384), _symm(0, frac(4, 5))),
                (w1, numpy.array([[lambda2, lambda2, 1 - 2*lambda2]])),
                (w2, numpy.array([[lambda3, lambda3, 1 - 2*lambda3]])),
                ]
        elif index == 9:
            self.degree = 3
            # ERR the article is missing the minus sign
            data = [
                (-frac(9, 32), _s3(symbolic)),
                (frac(25, 96), _s21(frac(1, 5))),
                ]
        else:
            assert index == 10
            self.degree = 3
            lambda1, lambda2 = [(16 + i*2*sqrt(14))/25 for i in [+1, -1]]
            w1, w2 = [(161 + i * 17*sqrt(14))/2688 for i in [+1, -1]]
            data = [
                (w2, _symm(lambda1, 0)),
                (w1, _symm(0, lambda2)),
                (
                    frac(25, 96),
                    numpy.array([[frac(2, 5), frac(2, 5), frac(1, 5)]])
                ),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        self.weights *= 2
        return


def _symm(a, b):
    c = 1 - a - b
    return numpy.array([
        [a, b, c],
        [b, a, c],
        ])
