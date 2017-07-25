# -*- coding: utf-8 -*-
#
from ..ncube import stroud


class StroudN(object):
    def __init__(self, index):
        self.name = 'Stroud(quad, {})'.format(index)
        w = stroud.Stroud(2, index)
        self.weights = w.weights
        self.points = w.points
        self.degree = w.degree
        return
