# -*- coding: utf-8 -*-
#
from .stroud import Stroud

from .. import helpers
from ..ncube import transform, integrate
from ..ncube import ncube_points as rectangle_points


def show(*args, **kwargs):
    import matplotlib.pyplot as plt

    plot(*args, **kwargs)
    plt.show()
    return


def plot(scheme, quad=rectangle_points([0.0, 1.0], [0.0, 1.0]), show_axes=False):
    """Shows the quadrature points on a given quad. The area of the disks
    around the points coincides with their weights.
    """
    import matplotlib.pyplot as plt

    plt.plot(quad[0][0], quad[1][0], "-k")
    plt.plot(quad[1][0], quad[1][1], "-k")
    plt.plot(quad[1][1], quad[0][1], "-k")
    plt.plot(quad[0][1], quad[0][0], "-k")

    plt.axis("equal")

    if not show_axes:
        plt.gca().set_axis_off()

    transformed_pts = transform(scheme.points.T, quad)

    vol = integrate(lambda x: 1.0, quad, Stroud("C2 1-1"))
    helpers.plot_disks(plt, transformed_pts, scheme.weights, vol)
    return
