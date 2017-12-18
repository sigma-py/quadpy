# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy

from .albrecht_collatz import AlbrechtCollatz
from .hammer_marlowe_stroud import HammerMarloweStroud

from ..helpers import untangle
from ..line_segment.gauss_legendre import GaussLegendre


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    def __init__(self, index):
        self.name = 'Stroud({})'.format(index)
        if index == 'T2 3-1':
            self.set_data(AlbrechtCollatz())
        elif index == 'T2 5-1':
            self.set_data(HammerMarloweStroud(5))
        else:
            # conical product Gauss
            assert index == 'T2 7-1'
            self.degree = 7

            gl4 = GaussLegendre(4)
            r = (gl4.points + 1) / 2
            A = gl4.weights / 2

            # Generate Gauss formula for int_0^1 (1-s) * f(s) ds.
            # ```
            # k = numpy.arange(8)
            # moments = 1 / (k**2 + 3*k + 2)
            # alpha, beta = orthopy.line.chebyshev(moments)
            # s, B = orthopy.line.schemes.custom(alpha, beta, mode='numpy')
            # ```
            s = numpy.array([
                5.710419611452533e-02,
                2.768430136381415e-01,
                5.835904323689318e-01,
                8.602401356562251e-01,
                ])
            B = numpy.array([
                1.355069134315012e-01,
                2.034645680102685e-01,
                1.298475476082247e-01,
                3.118097095000554e-02,
                ])

            data = [
                (2*A[i]*B[j], numpy.array([[
                    s[j], r[i]*(1-s[j]), (1-r[i])*(1-s[j])
                    ]]))
                for i in range(4)
                for j in range(4)
                ]

            self.bary, self.weights = untangle(data)
            self.points = self.bary[:, 1:]

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        self.bary = scheme.bary
        return
