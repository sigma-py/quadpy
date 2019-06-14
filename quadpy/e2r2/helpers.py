# -*- coding: utf-8 -*-
#
import collections
import numpy


E2r2Scheme = collections.namedtuple(
    "E2r2Scheme", ["name", "degree", "weights", "points"]
)


def _z():
    return numpy.array([[0.0, 0.0]])


def _s8(a, b):
    return numpy.array(
        [[+a, +b], [-a, +b], [+a, -b], [-a, -b], [+b, +a], [-b, +a], [+b, -a], [-b, -a]]
    )


def _s4(a):
    return numpy.array([[+a, +a], [-a, +a], [+a, -a], [-a, -a]])


def _s40(a):
    return numpy.array([[+a, 0.0], [-a, 0.0], [0.0, +a], [0.0, -a]])
