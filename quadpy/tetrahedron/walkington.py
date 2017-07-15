# -*- coding: utf-8 -*-
#
from ..simplex import walkington


class Walkington(object):
    def __init__(self, index):
        self.name = 'Walkington(tetrahedron, {})'.format(index)
        w = walkington.Walkington(3, index)
        self.weights = w.weights
        self.points = w.points
        self.degree = w.degree
        return
