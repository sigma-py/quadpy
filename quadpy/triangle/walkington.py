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
            self.weights = numpy.concatenate([
                numpy.full(1, fr(9, 40)),
                numpy.full(3, fr(31, 240) + sqrt(15) / 1200),
                numpy.full(3, fr(31, 240) - sqrt(15) / 1200),
                ])
            self.bary = numpy.concatenate([
                _c(),
                _xi1(fr(2, 7) + sqrt(15)/21),
                _xi1(fr(2, 7) - sqrt(15)/21),
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
