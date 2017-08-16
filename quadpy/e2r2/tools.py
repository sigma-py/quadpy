# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy

from .. import helpers


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(scheme, show_axes=True):
    ax = plt.gca()
    # change default range so that new disks will work
    plt.axis('equal')

    if not show_axes:
        ax.set_axis_off()

    helpers.plot_disks(plt, scheme.points, scheme.weights, numpy.pi)
    return


def integrate(f, rule, sumfun=helpers.kahan_sum):
    ff = numpy.array(f(rule.points.T))
    return sumfun(rule.weights * ff, axis=-1)
