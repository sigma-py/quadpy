# -*- coding: utf-8 -*-
#
'''
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
'''
from __future__ import division

from math import pi

from . import stroud_secrest

from ..helpers import untangle


_gen = {
    '5-1': stroud_secrest.vii,
    '5-2': stroud_secrest.viii,
    '5-3': stroud_secrest.ix,
    # '7-1': stroud_secrest.x,
    '7-2': stroud_secrest.xi_,
    }


class Stroud(object):
    keys = _gen.keys()

    def __init__(self, key):
        self.degree, data = _gen[key]()
        self.points, self.weights = untangle(data)
        self.weights *= 8 * pi
        return
