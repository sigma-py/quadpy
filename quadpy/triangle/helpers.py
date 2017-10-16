# -*- coding: utf-8 -*-
#
import numpy
from sympy import Rational


def _s3():
    return numpy.full((1, 3), Rational(1, 3))


def _s21(a):
    a = numpy.array(a)
    b = 1 - 2*a
    out = numpy.array([
        [a, a, b],
        [a, b, a],
        [b, a, a],
        ])
    return out


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
    c = 1.0 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        ])


def _collapse0(a):
    '''Collapse all dimensions of `a` except the first.
    '''
    return numpy.reshape(a, (a.shape[0], numpy.prod(a.shape[1:])))


def untangle2(data):
    bary = []
    weights = []

    if 's3' in data:
        data['s3'] = numpy.array(data['s3']).T
        bary.append(_s3().T)
        weights.append(numpy.tile(data['s3'][0], 1))

    if 's2' in data:
        data['s2'] = numpy.array(data['s2']).T
        s2_data = _s21(data['s2'][1])
        bary.append(_collapse0(s2_data))
        weights.append(numpy.tile(data['s2'][0], 3))

    if 's1' in data:
        data['s1'] = numpy.array(data['s1']).T
        s1_data = _s111ab(*data['s1'][1:])
        bary.append(_collapse0(s1_data))
        weights.append(numpy.tile(data['s1'][0], 6))

    bary = numpy.column_stack(bary).T
    weights = numpy.concatenate(weights)
    return bary, weights
