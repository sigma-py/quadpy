# -*- coding: utf-8 -*-
#
import numpy


def _s4():
    return numpy.array([[0.25, 0.25, 0.25, 0.25]])


def _s31(a):
    b = 1.0 - 3 * a
    return numpy.array([[a, a, a, b], [a, a, b, a], [a, b, a, a], [b, a, a, a]])


def _s22(a):
    b = 0.5 - a
    return numpy.array(
        [
            [a, a, b, b],
            [a, b, a, b],
            [b, a, a, b],
            [a, b, b, a],
            [b, a, b, a],
            [b, b, a, a],
        ]
    )


def _s211(a, b):
    c = 1.0 - 2 * a - b
    return numpy.array(
        [
            [a, a, b, c],
            [a, b, a, c],
            [b, a, a, c],
            [a, b, c, a],
            [b, a, c, a],
            [b, c, a, a],
            [a, a, c, b],
            [a, c, a, b],
            [c, a, a, b],
            [a, c, b, a],
            [c, a, b, a],
            [c, b, a, a],
        ]
    )


def _s1111(a, b, c):
    d = 1.0 - a - b - c
    return numpy.array(
        [
            [a, b, c, d],
            [a, b, d, c],
            [a, c, b, d],
            [a, c, d, b],
            [a, d, b, c],
            [a, d, c, b],
            [b, a, c, d],
            [b, a, d, c],
            [b, c, a, d],
            [b, c, d, a],
            [b, d, a, c],
            [b, d, c, a],
            [c, a, b, d],
            [c, a, d, b],
            [c, b, a, d],
            [c, b, d, a],
            [c, d, a, b],
            [c, d, b, a],
            [d, a, b, c],
            [d, a, c, b],
            [d, b, a, c],
            [d, b, c, a],
            [d, c, a, b],
            [d, c, b, a],
        ]
    )
