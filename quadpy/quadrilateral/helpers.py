# -*- coding: utf-8 -*-
#
import numpy


def _z():
    return numpy.array([[0, 0]])


def _symm_r_0(r):
    return numpy.array([[+r, 0], [-r, 0], [0, +r], [0, -r]])


def _symm_r0(r):
    z = numpy.zeros_like(r)
    return numpy.array([[+r, z], [-r, z], [z, +r], [z, -r]])


def _symm_s(s):
    return numpy.array([[+s, +s], [-s, +s], [+s, -s], [-s, -s]])


def _symm_s_t(s, t):
    return numpy.array(
        [[+s, +t], [-s, +t], [+s, -t], [-s, -t], [+t, +s], [-t, +s], [+t, -s], [-t, -s]]
    )


def _pm(s, t):
    return numpy.array([[+s, +t], [-s, -t]])


def _pm2(s, t):
    return numpy.array([[+s, +t], [-s, +t], [+s, -t], [-s, -t]])


def _pmx(x):
    z = numpy.zeros_like(x)
    return numpy.array([[+x, z], [-x, z]])


def _pmy(y):
    z = numpy.zeros_like(y)
    return numpy.array([[z, +y], [z, -y]])


def _collapse0(a):
    """Collapse all dimensions of `a` except the first.
    """
    return a.reshape(a.shape[0], -1)


def unroll(data, symbolic=False):
    bary = []
    weights = []

    if "zero" in data:
        d = numpy.array(data["zero"]).T
        bary.append(numpy.zeros((1, 2)))
        weights.append(numpy.tile(d[0], 1))

    if "symm_r0" in data:
        d = numpy.array(data["symm_r0"]).T
        r0_data = _symm_r0(d[1])
        r0_data = numpy.swapaxes(r0_data, 0, 1)
        bary.append(_collapse0(r0_data).T)
        weights.append(numpy.tile(d[0], 4))

    if "symm_s" in data:
        d = numpy.array(data["symm_s"]).T
        s_data = _symm_s(d[1])
        s_data = numpy.swapaxes(s_data, 0, 1)
        bary.append(_collapse0(s_data).T)
        weights.append(numpy.tile(d[0], 4))

    if "symm_s_t" in data:
        d = numpy.array(data["symm_s_t"]).T
        s_data = _symm_s_t(*d[1:])
        s_data = numpy.swapaxes(s_data, 0, 1)
        bary.append(_collapse0(s_data).T)
        weights.append(numpy.tile(d[0], 8))

    if "pm" in data:
        d = numpy.array(data["pm"]).T
        s_data = _pm(*d[1:])
        s_data = numpy.swapaxes(s_data, 0, 1)
        bary.append(_collapse0(s_data).T)
        weights.append(numpy.tile(d[0], 2))

    if "pm2" in data:
        d = numpy.array(data["pm2"]).T
        s_data = _pm2(*d[1:])
        s_data = numpy.swapaxes(s_data, 0, 1)
        bary.append(_collapse0(s_data).T)
        weights.append(numpy.tile(d[0], 4))

    if "pmx" in data:
        d = numpy.array(data["pmx"]).T
        s_data = _pmx(d[1])
        s_data = numpy.swapaxes(s_data, 0, 1)
        bary.append(_collapse0(s_data).T)
        weights.append(numpy.tile(d[0], 2))

    if "pmy" in data:
        d = numpy.array(data["pmy"]).T
        s_data = _pmy(d[1])
        s_data = numpy.swapaxes(s_data, 0, 1)
        bary.append(_collapse0(s_data).T)
        weights.append(numpy.tile(d[0], 2))

    bary = numpy.concatenate(bary)
    weights = numpy.concatenate(weights)
    return bary, weights
