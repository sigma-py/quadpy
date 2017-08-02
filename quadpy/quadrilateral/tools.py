# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers
from .stroud import Stroud
from ..ncube import transform, integrate


def rectangle_points(x, y):
    '''Given the end points of a rectangle aligned with the coordinate axes,
    this returns the corner points of the cube in the correct data structure.
    '''
    return numpy.moveaxis(
            numpy.array(numpy.meshgrid(x, y, indexing='ij')),
            0, -1
            )


def show(
        scheme,
        quad=rectangle_points([0.0, 1.0], [0.0, 1.0]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given quad. The area of the disks
    around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt

    plt.plot(quad[0][0], quad[1][0], '-k')
    plt.plot(quad[1][0], quad[1][1], '-k')
    plt.plot(quad[1][1], quad[0][1], '-k')
    plt.plot(quad[0][1], quad[0][0], '-k')

    plt.axis('equal')

    if not show_axes:
        plt.gca().set_axis_off()

    transformed_pts = transform(scheme.points.T, quad)

    vol = integrate(lambda x: 1.0, quad, Stroud(1))
    helpers.plot_disks(
        plt, transformed_pts, scheme.weights, vol
        )
    plt.show()
    return
