# -*- coding: utf-8 -*-
#
'''
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
'''
from __future__ import division

import numpy
import sympy

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

    def __init__(self, key, symbolic=False):
        self.name = 'Stround_E3r({})'.format(key)
        self.degree, data = _gen[key](symbolic)
        self.points, self.weights = untangle(data)
        pi = sympy.pi if symbolic else numpy.pi
        self.weights *= 8 * pi
        return
