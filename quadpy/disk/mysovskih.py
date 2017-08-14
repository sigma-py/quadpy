# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle, z


class Mysovskih(object):
    '''
    I.P. Mysovskih,
    On the construction of cubature formulas for the simplest regions,
    Z. Vychisl. Mat. i. Mat. Fiz. 4, 3-14, 1964.
    '''
    def __init__(self, alpha):
        self.degree = 4
        i = numpy.arange(5)
        b = numpy.sqrt((alpha + 4.0)/(alpha + 6.0))

        r = b * numpy.cos(0.4*i*numpy.pi)
        s = b * numpy.sin(0.4*i*numpy.pi)
        x = numpy.array([r, s]).T

        B0 = 4.0 / (alpha + 4)**2
        B1 = (alpha + 2.0)*(alpha + 6.0) / 5.0 / (alpha + 4.0)**2

        data = [
            (B0, z(2)),
            (B1, x),
            ]
        self.points, self.weights = untangle(data)
        self.weights *= numpy.pi
        return
