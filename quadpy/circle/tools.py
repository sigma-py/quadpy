# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy

from .. import helpers


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(scheme, show_axes=False):
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
    return


def integrate(f, center, radius, rule, sumfun=helpers.kahan_sum):
    flt = numpy.vectorize(float)
    center = numpy.array(center)
    rr = numpy.multiply.outer(radius, flt(rule.points))
    rr = numpy.swapaxes(rr, 0, -2)
    ff = numpy.array(f((rr + center).T))
    out = sumfun(flt(rule.weights) * ff, axis=-1)
    return radius * out
