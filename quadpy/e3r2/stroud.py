# -*- coding: utf-8 -*-
#
'''
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
'''
from __future__ import division

import numpy

from . import stroud_secrest

from ..sphere import stroud as sphere_stroud
from ..helpers import untangle


def _gen14_1(symbolic):
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
    # ```
    # n = 8
    # alpha = numpy.zeros(n)
    # beta = numpy.empty(n)
    # beta[0] = numpy.sqrt(numpy.pi)/2
    # beta[1::2] = numpy.arange(n//2) + 1.5
    # beta[2::2] = numpy.arange(n//2-1) + 1.0
    # points, weights = \
    #     orthopy.line.schemes.custom(alpha, beta, mode='numpy')
    # r = points[-4:]
    # A = weights[-4:]
    # ```
    r = numpy.array([
        7.235510187528402e-01,
        1.468553289216669e+00,
        2.266580584531844e+00,
        3.190993201781527e+00,
        ])
    A = numpy.array([
        2.265043732793035e-01,
        1.908084800858996e-01,
        2.539731378612040e-02,
        4.032955750550135e-04,
        ])

    spherical_scheme = sphere_stroud.Stroud('U3 14-1')
    v = spherical_scheme.points
    B = spherical_scheme.weights

    # Normalize the weights to 1
    B /= numpy.sqrt(numpy.pi) / 4

    data = [
        (A[i]*B[j], r[i] * numpy.array([v[j]]))
        for i in range(4)
        for j in range(72)
        ]

    return degree, data


# The boolean tells if the factor pi^{3/2} is already in the weights
_gen = {
    '5-1': stroud_secrest.vii,
    '5-2a': stroud_secrest.viiia,
    '5-2b': stroud_secrest.viiib,
    '5-3': stroud_secrest.ix,
    '7-1a': lambda symbolic: stroud_secrest.x(+1, symbolic),
    '7-1b': lambda symbolic: stroud_secrest.x(-1, symbolic),
    '7-2a': lambda symbolic: stroud_secrest.xi_(+1, symbolic),
    '7-2b': lambda symbolic: stroud_secrest.xi_(-1, symbolic),
    '14-1': _gen14_1
    }


class Stroud(object):
    keys = _gen.keys()

    def __init__(self, key, symbolic=False):
        self.name = 'Stround_E3r2({})'.format(key)
        self.degree, data = _gen[key](symbolic)
        self.points, self.weights = untangle(data)
        self.weights *= numpy.sqrt(numpy.pi)**3
        return
