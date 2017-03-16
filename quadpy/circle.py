# -*- coding: utf-8 -*-
#
import math
import numpy

from . import helpers


def show(scheme, show_axes=False):
    from matplotlib import pyplot as plt
    ax = plt.gca()
    # change default range so that new disks will work
    plt.axis('equal')
    ax.set_xlim((-1.5, 1.5))
    ax.set_ylim((-1.5, 1.5))

    if not show_axes:
        ax.set_axis_off()

    disk1 = plt.Circle((0, 0), 1, color='k', fill=False)
    ax.add_artist(disk1)

    # The total area is used to gauge the disk radii. This is only meaningful
    # for 2D manifolds, not for the circle. What we do instead is choose the
    # total_area such that the sum of the disk radii equals pi.
    total_area = numpy.pi**3 / len(scheme.weights)
    helpers.plot_disks(
        plt, scheme.points, scheme.weights, total_area
        )
    plt.show()
    return


def integrate(f, center, radius, rule, sum=helpers.kahan_sum):
    center = numpy.array(center)
    rr = numpy.multiply.outer(rule.points.T, radius)
    out = sum(
        (rule.weights * f(rr + center.T[:, None]).T).T,
        axis=0
        )
    return radius * out


class Equidistant(object):
    def __init__(self, n):
        self.weights = numpy.ones(n) * (2 * numpy.pi) / n
        self.points = numpy.column_stack([
            numpy.cos(2*numpy.pi * numpy.arange(n) / n),
            numpy.sin(2*numpy.pi * numpy.arange(n) / n),
            ])
        self.degree = n - 1
        return
