# -*- coding: utf-8 -*-
#
from .helpers import s3, TriangleScheme


def Centroid(symbolic=False):
    weights, bary = s3(1, symbolic)
    return TriangleScheme("Centroid rule", 1, weights, bary)
