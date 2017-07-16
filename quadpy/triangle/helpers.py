# -*- coding: utf-8 -*-
#
import numpy


def _s3():
    return numpy.array([
        [1.0/3.0, 1.0/3.0, 1.0/3.0]
        ])


def _s21(a):
    b = 1.0 - 2*a
    return numpy.array([
        [a, a, b],
        [a, b, a],
        [b, a, a],
        ])


def _s111(a, b):
    c = 1.0 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        [b, a, c],
        [c, b, a],
        [a, c, b],
        ])

