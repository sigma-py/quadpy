# -*- coding: utf-8 -*-
#
import numpy
from sympy import Rational


def _s3():
    return numpy.full((1, 3), Rational(1, 3))


def _s21(a):
    b = 1 - 2*a
    return numpy.array([
        [a, a, b],
        [a, b, a],
        [b, a, a],
        ])


def _s111(a, b):
    c = 1 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        [b, a, c],
        [c, b, a],
        [a, c, b],
        ])


def _rot(a, b):
    c = 1.0 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        ])
