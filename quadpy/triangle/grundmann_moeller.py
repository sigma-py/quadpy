# -*- coding: utf-8 -*-
#
from ..simplex import grundmann_moeller


class GrundmannMoeller(object):
    def __init__(self, s):
        self.name = 'GrundmannMöller(triangle, {})'.format(s)
        gm = grundmann_moeller.GrundmannMoeller(2, s)
        self.weights = gm.weights
        self.bary = gm.bary
        self.points = gm.points
        self.degree = gm.degree
        return
