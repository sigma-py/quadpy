# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy


def _s3(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x/y
    return numpy.full((1, 3), frac(1, 3))


def _s21(a):
    a = numpy.array(a)
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


def _s111ab(a, b):
    c = 1 - a - b
    out = numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        [b, a, c],
        [c, b, a],
        [a, c, b],
        ])
    out = numpy.swapaxes(out, 0, 1)
    return out


def _rot(a, b):
    c = 1 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        ])


def _rot_ab(a, b):
    c = 1 - a - b
    out = numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        ])
    out = numpy.swapaxes(out, 0, 1)
    return out


def _collapse0(a):
    '''Collapse all dimensions of `a` except the first.
    '''
    return numpy.reshape(a, (a.shape[0], numpy.prod(a.shape[1:])))


def untangle2(data, symbolic=False):
    bary = []
    weights = []

    if 's3' in data:
        d = numpy.array(data['s3']).T
        bary.append(_s3(symbolic).T)
        weights.append(numpy.tile(d[0], 1))

    if 's2' in data:
        d = numpy.array(data['s2']).T
        s2_data = _s21(d[1])
        bary.append(_collapse0(s2_data))
        weights.append(numpy.tile(d[0], 3))

    if 's1' in data:
        d = numpy.array(data['s1']).T
        s1_data = _s111ab(*d[1:])
        bary.append(_collapse0(s1_data))
        weights.append(numpy.tile(d[0], 6))

    if 'rot' in data:
        d = numpy.array(data['rot']).T
        rot_data = _rot_ab(*d[1:])
        bary.append(_collapse0(rot_data))
        weights.append(numpy.tile(d[0], 3))

    bary = numpy.column_stack(bary).T
    weights = numpy.concatenate(weights)
    return bary, weights
