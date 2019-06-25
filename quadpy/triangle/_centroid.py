# -*- coding: utf-8 -*-
#
from ._helpers import s3, TriangleScheme


def centroid(symbolic=False):
    weights, points = s3(1, symbolic)
    return TriangleScheme("Centroid rule", weights, points, 1)
