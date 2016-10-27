# -*- coding: utf-8 -*-
#
import math
import numpy

from . import helpers


def show(scheme):
    from matplotlib import pyplot as plt
    ax = plt.gca()
    # change default range so that new disks will work
    plt.axis('equal')
    ax.set_xlim((-1.5, 1.5))
    ax.set_ylim((-1.5, 1.5))

    disk1 = plt.Circle((0, 0), 1, color='k', fill=False)
    ax.add_artist(disk1)

    helpers.plot_disks(
        plt, scheme.points, scheme.weights, numpy.pi
        )
    return


def integrate(f, scheme):
    x = scheme.points.T
    return math.fsum(scheme.weights * f(x).T)


class Equidistant(object):
    def __init__(self, n):
        self.weights = numpy.ones(n) * (2 * numpy.pi) / n
        self.points = numpy.column_stack([
            numpy.cos(2*numpy.pi * numpy.arange(n) / n),
            numpy.sin(2*numpy.pi * numpy.arange(n) / n),
            ])
        self.degree = n - 1
        return
