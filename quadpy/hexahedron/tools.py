# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy

from .. import helpers
from ..ncube import transform, integrate
from ..ncube import ncube_points as cube_points


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(
        scheme,
        hexa=cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given hexahedron. The size of the
    balls around the points coincides with their weights.
    '''
    # pylint: disable=relative-import, unused-variable
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    if not show_axes:
        ax.set_axis_off()

    edges = _get_edges(hexa)

    for edge in edges:
        plt.plot(*edge, color='k', linestyle='-')

    transformed_pts = transform(scheme.points.T, hexa)

    vol = integrate(lambda x: 1.0, hexa, scheme)
    helpers.plot_spheres(
        plt, ax, transformed_pts, scheme.weights, vol
        )
    return


def show_mayavi(
        scheme,
        hexa=cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]),
        ):
    helpers.show_mayavi(
            transform(scheme.points.T, hexa),
            scheme.weights,
            integrate(lambda x: 1.0, hexa, scheme),
            _get_edges(hexa)
            )
    return


def show_vtk(
        scheme,
        hexa=cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]),
        ):
    helpers.show_vtk(
            transform(scheme.points.T, hexa),
            scheme.weights,
            integrate(lambda x: 1.0, hexa, scheme),
            _get_edges(hexa)
            )
    return


def _get_edges(hexa):
    edges = numpy.array([
        [hexa[0, 0, 0], hexa[1, 0, 0]],
        [hexa[1, 0, 0], hexa[1, 1, 0]],
        [hexa[1, 1, 0], hexa[0, 1, 0]],
        [hexa[0, 1, 0], hexa[0, 0, 0]],
        #
        [hexa[0, 0, 1], hexa[1, 0, 1]],
        [hexa[1, 0, 1], hexa[1, 1, 1]],
        [hexa[1, 1, 1], hexa[0, 1, 1]],
        [hexa[0, 1, 1], hexa[0, 0, 1]],
        #
        [hexa[0, 0, 0], hexa[0, 0, 1]],
        [hexa[1, 0, 0], hexa[1, 0, 1]],
        [hexa[1, 1, 0], hexa[1, 1, 1]],
        [hexa[0, 1, 0], hexa[0, 1, 1]],
        ])
    edges = numpy.moveaxis(edges, 1, 2)
    return edges
