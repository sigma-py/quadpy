# -*- coding: utf-8 -*-
#
import numpy

from ..simplex import walkington


class Walkington(object):
    def __init__(self, index):
        self.name = 'Walkington(triangle, {})'.format(index)
        if index == 'p5':
            self.degree = 5
            self.weights = numpy.concatenate([
                numpy.full(1, 9.0/40.0),
                numpy.full(3, 31.0/240.0 + numpy.sqrt(15) / 1200.0),
                numpy.full(3, 31.0/240.0 - numpy.sqrt(15) / 1200.0),
                ])
            self.bary = numpy.concatenate([
                _c(),
                _xi1(2.0/7.0 + numpy.sqrt(15.0)/21.0),
                _xi1(2.0/7.0 - numpy.sqrt(15.0)/21.0),
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
    return numpy.array([
        numpy.full(3, 1.0/3.0)
        ])


def _xi1(a):
    b = 1.0 - 2*a
    return numpy.array([
        [b, a, a],
        [a, b, a],
        [a, a, b],
        ])
