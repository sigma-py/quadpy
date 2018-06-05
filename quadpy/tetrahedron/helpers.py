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


def untangle2(data):
    print('untangle 2')
    print(data)
    print()

    bary = []
    weights = []

    if "s4" in data:
        assert len(data["s4"]) == 1
        w = numpy.array(data["s4"]).T
        bary.append(_s4())
        weights.append(w[0])

    if "s31" in data:
        d = numpy.array(data["s31"]).T
        s31_data = numpy.moveaxis(_s31(d[1]), 0, 1)
        bary.append(_collapse0(s31_data).T)
        weights.append(numpy.tile(d[0], 4))

    if "s22" in data:
        d = numpy.array(data["s22"]).T
        s22_data = numpy.moveaxis(_s22(d[1]), 0, 1)
        bary.append(_collapse0(s22_data).T)
        weights.append(numpy.tile(d[0], 6))

    if "s211" in data:
        d = numpy.array(data["s211"]).T
        s211_data = numpy.moveaxis(_s211(*d[1:]), 0, 1)
        bary.append(_collapse0(s211_data).T)
        weights.append(numpy.tile(d[0], 12))

    if "s1111" in data:
        # TODO
        assert False
        d = numpy.array(data["s1111"]).T
        s1111_data = _s1111(*d[1:])
        bary.append(_collapse0(s1111_data).T)
        weights.append(numpy.tile(d[0], 24))

    for b in bary:
        print(b.shape)

    bary = numpy.concatenate(bary)
    weights = numpy.concatenate(weights)

    return bary, weights


def _collapse0(a):
    """Collapse all dimensions of `a` except the first.
    """
    return numpy.reshape(a, (a.shape[0], numpy.prod(a.shape[1:])))
