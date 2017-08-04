# -*- coding: utf-8 -*-
#
import numpy

def z():
    return numpy.array([[0.0, 0.0, 0.0]])


def fs_r00(a):
    return numpy.array([
        [+a, 0.0, 0.0],
        [0.0, +a, 0.0],
        [0.0, 0.0, +a],
        [-a, 0.0, 0.0],
        [0.0, -a, 0.0],
        [0.0, 0.0, -a],
        ])


def fs_rr0(a):
    return numpy.array([
        [+a, +a, 0.0],
        [+a, 0.0, +a],
        [0.0, +a, +a],
        [+a, -a, 0.0],
        [+a, 0.0, -a],
        [0.0, +a, -a],
        [-a, +a, 0.0],
        [-a, 0.0, +a],
        [0.0, -a, +a],
        [-a, -a, 0.0],
        [-a, 0.0, -a],
        [0.0, -a, -a],
        ])


def pm_rrr(a):
    return numpy.array([
        [+a, +a, +a],
        [-a, +a, +a],
        [+a, -a, +a],
        [-a, -a, +a],
        [+a, +a, -a],
        [-a, +a, -a],
        [+a, -a, -a],
        [-a, -a, -a],
        ])
