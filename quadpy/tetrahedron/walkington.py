# -*- coding: utf-8 -*-
#
import numpy

from ..simplex import walkington


class Walkington(object):
    def __init__(self, index):
        self.name = 'Walkington(tetrahedron, {})'.format(index)

        if index == 'p5':
            self.degree = 5
            self.weights = 6 * numpy.concatenate([
                numpy.full(4, 0.018781320953002641800),
                numpy.full(4, 0.012248840519393658257),
                numpy.full(6, 0.0070910034628469110730),
                ])
            self.points = numpy.concatenate([
                _xi1(0.31088591926330060980),
                _xi1(0.092735250310891226402),
                _xi11(0.045503704125649649492),
                ])
            return

        # Default: scheme from general simplex
        w = walkington.Walkington(3, index)
        self.weights = w.weights
        self.points = w.points
        self.degree = w.degree
        return


def _xi1(a):
    b = 1.0 - 3*a
    return numpy.array([
        [b, a, a, a],
        [a, b, a, a],
        [a, a, b, a],
        [a, a, a, b],
        ])


def _xi11(a):
    b = (1.0 - 2*a) / 2.0
    return numpy.array([
        [b, b, a, a],
        [b, a, b, a],
        [b, a, a, b],
        [a, b, a, b],
        [a, a, b, b],
        [a, b, b, a],
        ])
