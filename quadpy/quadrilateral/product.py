# -*- coding: utf-8 -*-
#
import numpy


class Product(object):
    def __init__(self, scheme1d):
        self.schemes = \
            scheme1d if isinstance(scheme1d, list) \
            else 2 * [scheme1d]

        self.weights = numpy.outer(
            self.schemes[0].weights, self.schemes[1].weights
            ).flatten()
        self.points = numpy.dstack(numpy.meshgrid(
            self.schemes[0].points, self.schemes[1].points
            )).reshape(-1, 2)
        self.degree = min([s.degree for s in self.schemes])
        return
