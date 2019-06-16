# -*- coding: utf-8 -*-
#
from .helpers import s3, TriangleScheme


def Centroid(symbolic=False):
    weights, bary = s3(1, symbolic)
    points = bary[:, 1:]
    return TriangleScheme("Centroid rule", 1, weights, points, bary)
