# -*- coding: utf-8 -*-
#
'''
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
'''
import numpy
from sympy import sqrt, pi, Rational as fr, sin, cos

from . import rabinowitz_richter
from . import stroud_secrest
from ..helpers import untangle


def _gen4_1():
    pts = 2 * sqrt(5) * numpy.array([
        [cos(2*i*pi/5) for i in range(5)],
        [sin(2*i*pi/5) for i in range(5)],
        ]).T
    data = [
        (fr(7, 10), numpy.array([[0, 0]])),
        (fr(3, 50), pts),
        ]
    return 4, data


# The boolean tells whether the factor 2*pi is already in the weights
_gen = {
    '4-1': (_gen4_1, False),
    '5-1': (stroud_secrest.v, False),
    '7-1': (stroud_secrest.vi, False),
    '9-1': (rabinowitz_richter.gen1, True),
    '11-1': (rabinowitz_richter.gen2, True),
    '11-2': (rabinowitz_richter.gen3, True),
    # ERR misprint in Stroud copied from original article
    # '13-1': (rabinowitz_richter.gen4,
    '15-1': (rabinowitz_richter.gen5, True),
    }


class Stroud(object):
    keys = _gen.keys()

    def __init__(self, key):
        self.degree, data = _gen[key][0]()
        weights_contain_2pi = _gen[key][1]
        self.points, self.weights = untangle(data)
        if not weights_contain_2pi:
            self.weights *= 2 * pi
        return
