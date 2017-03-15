# -*- coding: utf-8 -*-
#
import math
import numpy


def plot_disks(plt, pts, weights, total_area):
    '''Plot a circles at quadrature points according to weights.
    '''
    sum_weights = math.fsum(weights)
    for tp, weight in zip(pts, weights):
        # use matplotlib 2.0's color scheme
        color = '1f77b4' if weight >= 0 else 'd62728'
        # highlight circle center
        plt.plot([tp[0]], [tp[1]], '.' + color)
        # Choose radius such that the sum of areas of the circles equals
        # total_area.
        radius = math.sqrt(abs(weight)/sum_weights * total_area/math.pi)
        circ = plt.Circle((tp[0], tp[1]), radius, color=color, alpha=0.5)
        plt.gca().add_artist(circ)

    a = 1.3
    plt.gca().set_xlim(-a, +a)
    plt.gca().set_ylim(-a, +a)

    return


def plot_balls(
        plt, ax, pts, weights, total_volume,
        xmin, xmax, ymin, ymax, zmin, zmax
        ):
    phi, theta = numpy.mgrid[0:numpy.pi:101j, 0:2*numpy.pi:101j]
    x = numpy.sin(phi)*numpy.cos(theta)
    y = numpy.sin(phi)*numpy.sin(theta)
    z = numpy.cos(phi)

    alpha = 0.3

    sum_weights = math.fsum(weights)
    for tp, weight in zip(pts, weights):
        color = '1f77b4' if weight >= 0 else 'd62728'
        # highlight ball center
        plt.plot([tp[0]], [tp[1]], [tp[2]], '.' + color)

        # Choose radius such that the sum of volumes of the balls equals
        # total_volume.
        r = (
            abs(weight)/sum_weights * total_volume/(4.0/3.0 * numpy.pi)
            )**(1.0/3.0)

        ax.plot_surface(
            r*x + tp[0], r*y + tp[1], r*z + tp[2],
            color=color,
            alpha=alpha,
            linewidth=0
            )

    # http://stackoverflow.com/a/21765085/353337
    alpha = 1.3
    max_range = alpha * 0.5 * numpy.array([
        xmax - xmin,
        ymax - ymin,
        zmax - zmin,
        ]).max()
    mid_x = 0.5 * (xmax + xmin)
    mid_y = 0.5 * (ymax + ymin)
    mid_z = 0.5 * (zmax + zmin)
    #
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_axis_off()

    return
