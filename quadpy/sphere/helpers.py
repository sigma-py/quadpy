# -*- coding: utf-8 -*-
#
import numpy
import sympy


def cartesian_to_spherical(X):
    return numpy.stack([
        numpy.arctan2(X[:, 1], X[:, 0]),
        numpy.arccos(X[:, 2])
        ], axis=1)


def _atan2_0(X):
    '''Like sympy.atan2, but return 0 for x=y=0. Mathematically, the value is
    undefined, so sympy returns NaN, but for the sake of the coordinate
    conversion, its value doesn't matter. NaNs, however, produce NaNs down the
    line.
    '''
    out = numpy.array([sympy.atan2(X[k, 1], X[k, 0]) for k in range(len(X))])
    out[out == sympy.nan] = 0
    return out


def cartesian_to_spherical_sympy(X):
    vacos = numpy.vectorize(sympy.acos)
    return numpy.stack([
        _atan2_0(X),
        vacos(X[:, 2])
        ], axis=1)
