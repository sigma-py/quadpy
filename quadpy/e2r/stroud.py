# -*- coding: utf-8 -*-
#
'''
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
'''
from __future__ import division

from math import sqrt, pi, fsum

import numpy

from . import rabinowitz_richter
from . import stroud_secrest

from ..helpers import untangle


def _gen4_1():
    i = numpy.arange(5)
    pts = 2 * sqrt(5) * numpy.array([
        numpy.cos(2*i*pi/5),
        numpy.sin(2*i*pi/5),
        ]).T
    data = [
        (7/10, numpy.array([[0.0, 0.0]])),
        (3/50, pts),
        ]
    return 4, data


_gen = {
    '4-1': _gen4_1,
    '5-1': stroud_secrest.v,
    '7-1': stroud_secrest.vi,
    '9-1': rabinowitz_richter.gen1,
    '11-1': rabinowitz_richter.gen2,
    '11-2': rabinowitz_richter.gen3,
    # ERR misprint in Stroud copied from original article
    # '13-1': rabinowitz_richter.gen4,
    '15-1': rabinowitz_richter.gen5,
    }


class Stroud(object):
    keys = _gen.keys()

    def __init__(self, key):
        self.degree, data = _gen[key]()
        self.points, self.weights = untangle(data)
        self.weights /= fsum(self.weights)
        self.weights *= 2 * pi
        return
