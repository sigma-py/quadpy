# -*- coding: utf-8 -*-
#

def _z():
    return numpy.array([0.0, 0.0])


def _symm_r_0(r):
    return numpy.array([
        [+r, 0.0],
        [-r, 0.0],
        [0.0, +r],
        [0.0, -r],
        ])


def _symm_s(s):
    return numpy.array([
        [+s, +s],
        [-s, +s],
        [+s, -s],
        [-s, -s],
        ])


def _symm_s_t(s, t):
    return numpy.array([
        [+s, +t],
        [-s, +t],
        [+s, -t],
        [-s, -t],
        [+t, +s],
        [-t, +s],
        [+t, -s],
        [-t, -s],
        ]_
