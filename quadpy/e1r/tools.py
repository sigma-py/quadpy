# -*- coding: utf-8 -*-
#
from math import fsum

import matplotlib.pyplot as plt
import numpy

from .. import helpers


def integrate(f, rule, sumfun=helpers.kahan_sum):
    ff = numpy.array(f(rule.points.T))
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
    helpers.plot_disks_1d(plt, pts, scheme.weights, total_area=1.0)
    return
