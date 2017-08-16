# -*- coding: utf-8 -*-
#
import math

import matplotlib.pyplot as plt
import numpy

from .. import helpers


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(scheme, show_axes=True):
    ax = plt.gca()
    plt.axis('equal')

    if not show_axes:
        ax.set_axis_off()

    n = 2
    I0 = 2*math.factorial(n-1)*math.pi**(0.5*n) / math.gamma(0.5*n)

    helpers.plot_disks(
        plt, scheme.points, scheme.weights, I0
        )
    return


def integrate(f, rule, sumfun=helpers.kahan_sum):
    ff = numpy.array(f(rule.points.T))
    return sumfun(rule.weights * ff, axis=-1)
