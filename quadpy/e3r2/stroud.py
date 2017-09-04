# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import orthopy

from .stroud_secrest import StroudSecrest

from ..sphere import stroud as sphere_stroud
from ..helpers import untangle


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index):
        if index == 'E3r2 5-1':
            self.set_data(StroudSecrest('VII'))
        elif index == 'E3r2 5-2a':
            self.set_data(StroudSecrest('VIIIa'))
        elif index == 'E3r2 5-2b':
            self.set_data(StroudSecrest('VIIIb'))
        elif index == 'E3r2 5-3':
            self.set_data(StroudSecrest('IX'))
        elif index == 'E3r2 7-1a':
            self.set_data(StroudSecrest('Xa'))
        elif index == 'E3r2 7-1b':
            self.set_data(StroudSecrest('Xb'))
        elif index == 'E3r2 7-2a':
            self.set_data(StroudSecrest('XIa'))
        elif index == 'E3r2 7-2b':
            self.set_data(StroudSecrest('XIb'))
        else:
            assert index == 'E3r2 14-1'
            self.degree = 14
            # Get the moments corresponding to monomials and the weight
            # function omega(x) = x^2 * exp(-x^2):
            #
            #    int_{-infty}^{infty} x^2 exp(-x^2) x^k dx \
            #
            #           / 0 for k odd,
            #        = {
            #           \ Gamma((k+3)/2) if k even
            #
            # In this particular case, we don't need to compute the recurrence
            # coefficients numerically, but they are given analytically.
            n = 8
            alpha = numpy.zeros(n)
            beta = numpy.empty(n)
            beta[0] = numpy.sqrt(numpy.pi)/2
            beta[1::2] = numpy.arange(n//2) + 1.5
            beta[2::2] = numpy.arange(n//2-1) + 1.0

            points, weights = orthopy.schemes.custom(alpha, beta, mode='numpy')

            r = points[-4:]
            A = weights[-4:]

            spherical_scheme = sphere_stroud.Stroud('U3 14-1')
            v = spherical_scheme.points
            B = spherical_scheme.weights

            data = [
                (A[i]*B[j], r[i] * numpy.array([v[j]]))
                for i in range(4)
                for j in range(72)
                ]

            self.points, self.weights = untangle(data)
            self.weights *= 4 * numpy.pi

        return

    def set_data(self, scheme):
        self.degree = scheme.degree
        self.weights = scheme.weights
        self.points = scheme.points
        return
