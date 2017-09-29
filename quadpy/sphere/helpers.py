# -*- coding: utf-8 -*-
#
import numpy
import sympy


def cartesian_to_spherical(X):
    return numpy.stack([
        numpy.arctan2(X[:, 1], X[:, 0]),
        numpy.arccos(X[:, 2])
        ], axis=1)


def cartesian_to_spherical_sympy(X):
    vacos = numpy.vectorize(sympy.acos)
    return numpy.stack([
        [sympy.atan2(X[k, 1], X[k, 0]) for k in range(len(X))],
        vacos(X[:, 2])
        ], axis=1)
