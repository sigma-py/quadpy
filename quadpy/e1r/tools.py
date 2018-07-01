# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def integrate(f, rule, dot=numpy.dot):
    return dot(f(rule.points.T), rule.weights)


def show(*args, **kwargs):
    import matplotlib.pyplot as plt
    plot(*args, **kwargs)
    plt.show()
    return


def plot(scheme):
    import matplotlib.pyplot as plt
    plt.axis("equal")
    m = 1.1 * numpy.max(scheme.points)
    plt.plot([0, m], [0, 0], color="k")
    pts = numpy.column_stack([scheme.points, numpy.zeros(len(scheme.points))])
    helpers.plot_disks_1d(plt, pts, scheme.weights, total_area=1.0)
    return
