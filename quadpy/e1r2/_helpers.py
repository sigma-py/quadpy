# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import plot_disks_1d


class E1r2Scheme(object):
    def __init__(self, name, weights, points, degree):
        self.name = name
        self.weights = weights
        self.points = points
        self.degree = degree
        return

    def integrate(self, f, dot=numpy.dot):
        return dot(f(self.points.T), self.weights)

    def show(self, *args, **kwargs):
        import matplotlib.pyplot as plt

        self.plot(*args, **kwargs)
        plt.show()
        return

    def plot(self):
        import matplotlib.pyplot as plt

        plt.axis("equal")
        m = 1.1 * numpy.max(self.points)
        plt.plot([-m, +m], [0, 0], color="k")
        pts = numpy.column_stack([self.points, numpy.zeros(len(self.points))])
        plot_disks_1d(plt, pts, self.weights, total_area=1.0)
        return
