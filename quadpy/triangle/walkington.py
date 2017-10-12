# -*- coding: utf-8 -*-
#
import numpy
from sympy import Rational as fr, sqrt

from ..nsimplex import walkington


class Walkington(object):
    def __init__(self, index):
        self.name = 'Walkington(triangle, {})'.format(index)
        if index == 'p5':
            self.degree = 5

            a1, a2 = [(155 + i * sqrt(15))/1200 for i in [+1, -1]]
            self.weights = numpy.concatenate([
                numpy.full(1, fr(9, 40)),
                numpy.full(3, a1),
                numpy.full(3, a2),
                ])

            x1, x2 = [(6 + i * sqrt(15))/21 for i in [+1, -1]]
            self.bary = numpy.concatenate([
                _c(),
                _xi1(x1),
                _xi1(x2),
                ])
            self.points = self.bary[:, 1:]
            return

        # Default: scheme from general simplex
        w = walkington.Walkington(2, index)
        self.weights = w.weights
        self.bary = w.bary
        self.points = w.points
        self.degree = w.degree
        return


def _c():
    return numpy.full((1, 3), fr(1, 3))


def _xi1(a):
    b = 1 - 2*a
    return numpy.array([
        [b, a, a],
        [a, b, a],
        [a, a, b],
        ])
