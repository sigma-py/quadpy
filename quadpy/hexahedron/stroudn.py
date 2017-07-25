# -*- coding: utf-8 -*-
#
from ..ncube import stroud


class StroudN(object):
    def __init__(self, index):
        self.name = 'Stroud(hex, {})'.format(index)
        w = stroud.Stroud(3, index)
        self.weights = w.weights
        self.points = w.points
        self.degree = w.degree
        return
