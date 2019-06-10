# -*- coding: utf-8 -*-
#
import numpy


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


def _pm(a, b):
    return numpy.array([[+a, +b], [-a, +b], [+a, -b], [-a, -b]])


def _pmx(x):
    return numpy.array([[+x, 0], [-x, 0]])


def _pmy(y):
    return numpy.array([[0, +y], [0, -y]])
