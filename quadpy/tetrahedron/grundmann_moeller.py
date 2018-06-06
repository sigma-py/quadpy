# -*- coding: utf-8 -*-
#
from ..nsimplex import grundmann_moeller


class GrundmannMoeller(object):
    def __init__(self, s):
        self.name = "GrundmannMöller(tetrahedron, {})".format(s)
        gm = grundmann_moeller.GrundmannMoeller(3, s)
        self.weights = gm.weights
        self.points = gm.points
        self.degree = gm.degree
        return
