# -*- coding: utf-8 -*-
#
import numpy
from .helpers import _s3, TriangleScheme


def Centroid(symbolic=False):
    weights = numpy.array([1])
    bary = _s3(symbolic)
    points = bary[:, 1:]
    return TriangleScheme("Centroid rule", 1, weights, points, bary)
