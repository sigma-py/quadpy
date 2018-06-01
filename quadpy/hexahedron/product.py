# -*- coding: utf-8 -*-
#
import numpy


class Product(object):
    def __init__(self, scheme1d):
        self.schemes = scheme1d if isinstance(scheme1d, list) else 3 * [scheme1d]

        wy, wz, wx = numpy.meshgrid(
            self.schemes[0].weights, self.schemes[1].weights, self.schemes[2].weights
        )
        weights = numpy.vstack([wx.flatten(), wy.flatten(), wz.flatten()]).T
        self.weights = numpy.prod(weights, axis=1)
        # the order, yeah...
        y, z, x = numpy.meshgrid(
            self.schemes[0].points, self.schemes[1].points, self.schemes[2].points
        )
        self.points = numpy.vstack([x.flatten(), y.flatten(), z.flatten()]).T

        self.degree = min([s.degree for s in self.schemes])
        return
