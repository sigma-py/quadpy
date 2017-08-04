# -*- coding: utf-8 -*-
#
import numpy


def fs00(a):
    return numpy.array([
        [+a, 0.0, 0.0],
        [0.0, +a, 0.0],
        [0.0, 0.0, +a],
        [-a, 0.0, 0.0],
        [0.0, -a, 0.0],
        [0.0, 0.0, -a],
        ])
