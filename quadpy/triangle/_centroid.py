# -*- coding: utf-8 -*-
#
from ._helpers import TriangleScheme, s3


def centroid(symbolic=False):
    weights, points = s3(1, symbolic)
    return TriangleScheme("Centroid rule", weights, points, 1)
