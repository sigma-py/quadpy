# -*- coding: utf-8 -*-
#
"""
Arthur Stroud,
Approximate Calculation of Multiple Integrals,
Prentice Hall, 1971.
"""
from __future__ import division

import numpy

from .stroud_secrest import StroudSecrest

from ..sphere import _stroud as sphere_stroud
from ..helpers import untangle
from .helpers import E3r2Scheme


def stroud_14_1(symbolic=False):
    # Get the moments corresponding to monomials and the weight function omega(x) = x^2
    # * exp(-x^2):
    #
    #    int_{-infty}^{infty} x^2 exp(-x^2) x^k dx \
    #
    #           / 0 for k odd,
    #        = {
    #           \ Gamma((k+3)/2) if k even
    #
    # In this particular case, we don't need to compute the recurrence coefficients
    # numerically, but they are given analytically.
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
    r = numpy.array(
        [
            7.235510187528402e-01,
            1.468553289216669e00,
            2.266580584531844e00,
            3.190993201781527e00,
        ]
    )
    A = numpy.array(
        [
            2.265043732793035e-01,
            1.908084800858996e-01,
            2.539731378612040e-02,
            4.032955750550135e-04,
        ]
    )

    spherical_scheme = sphere_stroud.Stroud("U3 14-1")
    v = spherical_scheme.points
    B = spherical_scheme.weights

    # Normalize the weights to 1
    B /= numpy.sqrt(numpy.pi) / 4

    data = [
        (A[i] * B[j], r[i] * numpy.array([v[j]])) for i in range(4) for j in range(72)
    ]

    points, weights = untangle(data)
    weights *= numpy.sqrt(numpy.pi) ** 3
    return E3r2Scheme("Stroud 14-1", 14, weights, points)


Stroud = {
    "5-1": StroudSecrest["VII"],
    "5-2a": StroudSecrest["VIIIa"],
    "5-2b": StroudSecrest["VIIIb"],
    "5-3": StroudSecrest["IX"],
    "7-1a": StroudSecrest["Xa"],
    "7-1b": StroudSecrest["Xb"],
    "7-2a": StroudSecrest["XIa"],
    "7-2b": StroudSecrest["XIb"],
    "14-1": stroud_14_1,
}
