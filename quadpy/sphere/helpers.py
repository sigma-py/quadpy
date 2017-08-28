# -*- coding: utf-8 -*-
#
import numpy


def cartesian_to_spherical(X):
    return numpy.stack([
        numpy.arctan2(X[:, 1], X[:, 0]),
        numpy.arccos(X[:, 2])
        ], axis=1)
