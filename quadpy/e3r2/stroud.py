# -*- coding: utf-8 -*-
#
'''
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
'''
from __future__ import division

from math import fsum

import numpy
import orthopy

from . import stroud_secrest

from ..sphere import stroud as sphere_stroud
from ..helpers import untangle


def _gen14_1():
    degree = 14
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

    return degree, data


_gen = {
    '5-1': stroud_secrest.vii,
    '5-2a': stroud_secrest.viiia,
    '5-2b': stroud_secrest.viiib,
    '5-3': stroud_secrest.ix,
    '7-1a': lambda: stroud_secrest.x(+1),
    '7-1b': lambda: stroud_secrest.x(-1),
    '7-2a': lambda: stroud_secrest.xi_(+1),
    '7-2b': lambda: stroud_secrest.xi_(-1),
    '14-1': _gen14_1,
    }


class Stroud(object):
    keys = _gen.keys()

    def __init__(self, key):
        self.degree, data = _gen[key]()
        self.points, self.weights = untangle(data)
        self.weights /= fsum(self.weights)
        self.weights *= numpy.sqrt(numpy.pi)**3
        return
