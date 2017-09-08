# -*- coding: utf-8 -*-
#
from math import fsum

import matplotlib.pyplot as plt
import numpy

from .. import helpers


def integrate(f, rule, sumfun=helpers.kahan_sum):
    ff = numpy.array(f(rule.points.T))
    print(ff.shape)
    print(rule.points.shape)
    print(rule.weights.shape)
    return sumfun(rule.weights * ff, axis=-1)


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(scheme):
    plt.axis('equal')

    m = 1.1 * numpy.max(scheme.points)
    plt.plot([0, m], [0, 0], color='k')

    pts = numpy.column_stack([scheme.points, numpy.zeros(len(scheme.points))])

    # The total area is used to gauge the disk radii. This is only meaningful
    # for 2D manifolds, not for the circle. What we do instead is choose the
    # total_area such that the sum of the disk radii equals b-a.
    length = m
    total_area = 0.25 * length**2 * numpy.pi * fsum(scheme.weights) \
        / fsum(numpy.sqrt(abs(scheme.weights)))**2

    helpers.plot_disks(
        plt, pts, scheme.weights, total_area
        )
    return
