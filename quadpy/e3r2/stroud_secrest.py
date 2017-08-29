# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from ..helpers import untangle, pm_roll


class StroudSecrest(object):
    '''
    A.H. Stroud and D. Secrest,
    Approximate integration formulas for certain spherically symmetric regions,
    Math. Comp. 17 (1963), 105-135,
    <https://doi.org/10.1090/S0025-5718-1963-0161473-0>.
    '''
    def __init__(self, index):
        if index == 'VII':
            self.degree = 5

            pm = numpy.array([+1, -1])

            # article:
            # nu, xi = numpy.sqrt((15 + pm * 3*numpy.sqrt(5)))
            # A = 3/5
            # B = 1/30

            # book:
            nu, xi = numpy.sqrt((5 - pm * numpy.sqrt(5)) / 4)
            A = 2/5
            B = 1/20

            data = [
                (A, numpy.array([[0.0, 0.0, 0.0]])),
                (B, pm_roll(3, [nu, xi])),
                ]
        else:
            assert False

        self.points, self.weights = untangle(data)
        self.weights *= numpy.sqrt(numpy.pi)**3
        return
